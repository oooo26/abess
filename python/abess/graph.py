import numbers
import numpy as np
from scipy.sparse import coo_matrix
from sklearn.utils.validation import check_array
from .cabess import *
from .bess_base import bess_base


def fix_docs(cls):
    # inherit the document from base class
    index = cls.__doc__.find("Examples\n    --------\n")
    if index != -1:
        cls.__doc__ = cls.__doc__[:index] + \
            cls.__bases__[0].__doc__ + cls.__doc__[index:]

    # for name, func in vars(cls).items():
    #     if isinstance(func, types.FunctionType):
    #         # print(str(func) +  'needs doc')
    #         for parent in cls.__bases__:
    #             parfunc = getattr(parent, name, None)
    #             if parfunc and getattr(parfunc, '__doc__', None):
    #                 func.__doc__ = parfunc.__doc__ + func.__doc__
    return cls


@fix_docs
class abessDAG(bess_base):
    """
    Adaptive Best-Subset Selection(ABESS) algorithm for principal component analysis.

    Parameters
    ----------
    splicing_type: {0, 1}, optional
        The type of splicing in `fit()` (in Algorithm.h).
        "0" for decreasing by half, "1" for decresing by one.
        Default: splicing_type = 1.

    Examples
    --------
    >>> ### Sparsity known
    >>>
    >>> from abess.pca import abessPCA
    >>> import numpy as np
    >>> np.random.seed(12345)
    >>> model = abessPCA(support_size = 10)
    >>>
    >>> ### X known
    >>> X = np.random.randn(100, 50)
    >>> model.fit(X)
    >>> print(model.coef_)
    >>>
    >>> ### X unknown, but Sigma known
    >>> model.fit(Sigma = np.cov(X.T))
    >>> print(model.coef_)
    """

    def __init__(self, max_iter=20, exchange_num=5, path_type="seq", is_warm_start=True, support_size=None, s_min=None, s_max=None,
                 ic_type="ebic", ic_coef=1.0, cv=1, screening_size=-1,
                 always_select=None,
                 thread=1,
                 sparse_matrix=False,
                 splicing_type=1
                 ):
        super().__init__(
            algorithm_type="abess", model_type="DAG", normalize_type=1, path_type=path_type, max_iter=max_iter, exchange_num=exchange_num,
            is_warm_start=is_warm_start, support_size=support_size, s_min=s_min, s_max=s_max,
            ic_type=ic_type, ic_coef=ic_coef, cv=cv, screening_size=screening_size,
            always_select=always_select,
            thread=thread,
            sparse_matrix=sparse_matrix,
            splicing_type=splicing_type
        )

    def fit(self, X=None, is_normal=False, A_init=None):
        """
        The fit function is used to transfer the information of data and return the fit result.

        Parameters
        ----------
        X : array-like of shape (n_samples, p_features)
            Training data
        is_normal : bool, optional
            whether normalize the variables array before fitting the algorithm.
            Default: is_normal=False.
        weight : array-like of shape (n_samples,)
            Individual weights for each sample. Only used for is_weight=True.
            Default is 1 for each observation.
        group : int, optional
            The group index for each variable.
            Default: group = \\code{numpy.ones(p)}.
        Sigma : array-like of shape (n_features, n_features), optional
            Sample covariance matrix.
            For PCA, it can be given as input, instead of X. But if X is given, Sigma will be set to \\code{np.cov(X.T)}.
            Default: Sigma = \\code{np.cov(X.T)}.
        number : int, optional
            Indicates the number of PCs returned.
            Default: 1
        n : int, optional
            Sample size. If X is given, it would be X.shape[0]; if Sigma is given, it would be 1 by default.
            Default: X.shape[0] or 1.
        """

        # Input check
        if isinstance(X, (list, np.ndarray, np.matrix, coo_matrix)):
            if isinstance(X, coo_matrix):
                self.sparse_matrix = True
            X = check_array(X, accept_sparse=True)

            n = X.shape[0]
            p = X.shape[1]
            X = X - X.mean(axis=0)
            Sigma = np.cov(X.T)
            self.n_features_in_ = p
        else:
            raise ValueError("X should be given.")

        # # Algorithm_type
        # if self.algorithm_type == "abess":
        #     algorithm_type_int = 6
        # else:
        #     raise ValueError("algorithm_type should not be " +
        #                      str(self.algorithm_type))

        # for DAG,
        # model_type_int = 11
        path_type_int = 1

        # Ic_type
        if self.ic_type == "aic":
            ic_type_int = 1
        elif self.ic_type == "bic":
            ic_type_int = 2
        elif self.ic_type == "gic":
            ic_type_int = 3
        elif self.ic_type == "ebic":
            ic_type_int = 4
        else:
            raise ValueError(
                "ic_type should be \"aic\", \"bic\", \"ebic\" or \"gic\"")

        # cv
        if (not isinstance(self.cv, int) or self.cv <= 0):
            raise ValueError("cv should be an positive integer.")
        if self.cv > n:
            raise ValueError("cv should be smaller than n.")

        # Group
        g_index = list(range(p*p))
        # if group is None:
        #     g_index = list(range(p))
        # else:
        #     group = np.array(group)
        #     if group.ndim > 1:
        #         raise ValueError("group should be an 1D array of integers.")
        #     if group.size != p:
        #         raise ValueError(
        #             "The length of group should be equal to X.shape[1].")
        #     g_index = []
        #     group.sort()
        #     group_set = list(set(group))
        #     j = 0
        #     for i in group_set:
        #         while group[j] != i:
        #             j += 1
        #         g_index.append(j)

        # path parameter (note that: path_type_int = 1)
        if self.support_size is None:
            support_sizes = np.arange(0, p)
        else:
            if isinstance(self.support_size, (numbers.Integral)):
                support_sizes = [self.support_size]
            elif (len(self.support_size.shape) != 1):
                raise ValueError(
                    "`support_size` should be 1-dimension")
            elif self.support_size.max() > p*(p-1)/2:
                raise ValueError(
                    "`support_size` should not larger than p*(p-1)/2")
            else:
                support_sizes = self.support_size
        support_sizes = np.array(support_sizes).astype('int32')

        # screening
        if self.screening_size != -1:
            if self.screening_size == 0:
                self.screening_size = min(
                    p, int(n / (np.log(np.log(n)) * np.log(p))))
            elif self.screening_size > p:
                raise ValueError(
                    "screening size should be smaller than X.shape[1].")
            elif self.screening_size < max(support_sizes):
                raise ValueError(
                    "screening size should be more than max(support_size).")

        # unused
        new_s_min = 0
        new_s_max = 0
        cv_fold_id = np.array([], dtype="int32")

        # Exchange_num
        if (not isinstance(self.exchange_num, int) or self.exchange_num <= 0):
            raise ValueError("exchange_num should be an positive integer.")

        # Thread
        if (not isinstance(self.thread, int) or self.thread < 0):
            raise ValueError(
                "thread should be positive number or 0 (maximum supported by your device).")

        # Splicing type
        if self.splicing_type not in (0, 1):
            raise ValueError("splicing type should be 0 or 1.")

        # Important_search
        if (not isinstance(self.important_search, int)
                or self.important_search < 0):
            raise ValueError(
                "important_search should be a non-negative number.")

        # A_init
        if A_init is None:
            A_init = np.array([], dtype="int32")
        else:
            A_init = np.array(A_init, dtype="int32")
            if A_init.ndim > 1:
                raise ValueError(
                    "The initial active set should be an 1D array of integers.")
            if (A_init.min() < 0 or A_init.max() >= p*p):
                raise ValueError(
                    "A_init contains wrong index.")

        # Sparse X
        if self.sparse_matrix:
            if not isinstance(X, type(coo_matrix((1, 1)))):
                # print("sparse matrix 1")
                nonzero = 0
                tmp = np.zeros([X.shape[0] * X.shape[1], 3])
                for j in range(X.shape[1]):
                    for i in range(X.shape[0]):
                        if X[i, j] != 0.:
                            tmp[nonzero, :] = np.array([X[i, j], i, j])
                            nonzero += 1
                X = tmp[:nonzero, :]
            else:
                # print("sparse matrix 2")
                tmp = np.zeros([len(X.data), 3])
                tmp[:, 1] = X.row
                tmp[:, 2] = X.col
                tmp[:, 0] = X.data

                ind = np.lexsort((tmp[:, 2], tmp[:, 1]))
                X = tmp[ind, :]

        # normalize
        normalize = 0
        if is_normal:
            normalize = self.normalize_type

        # always_select
        if self.always_select is None:
            self.always_select = []

        # wrap with cpp
        # weight = np.ones(n)
        result = pywrap_DAG(X, 
                            n, p, normalize, 
                            self.max_iter, self.exchange_num,
                            path_type_int, self.is_warm_start,
                            ic_type_int, self.ic_coef, self.cv,
                            g_index,
                            support_sizes,
                            cv_fold_id,
                            new_s_min, new_s_max,
                            self.screening_size,
                            self.always_select,
                            self.early_stop,
                            self.thread,
                            self.sparse_matrix,
                            self.splicing_type,
                            self.important_search,
                            A_init,
                            p, 1,
                            1, 1
                            )

        self.coef_ = result[0]
        return self
