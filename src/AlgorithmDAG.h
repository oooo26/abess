#ifndef SRC_ALGORITHMDAG_H
#define SRC_ALGORITHMDAG_H

#include <Spectra/SymEigsSolver.h>

#include <cmath>

#include "Algorithm.h"

using namespace Spectra;

template <class T4>
class abessDAG : public Algorithm<Eigen::VectorXd, Eigen::VectorXd, double, T4> {
   public:
    // MatrixXd Adj;  // adjacency matrix
    // MatrixXd Sigma;

    abessDAG(int algorithm_type, int model_type, int max_iter = 30, int primary_model_fit_max_iter = 10,
             double primary_model_fit_epsilon = 1e-8, bool warm_start = true, int exchange_num = 5,
             Eigen::VectorXi always_select = Eigen::VectorXi::Zero(0), int splicing_type = 1, int sub_search = 0)
        : Algorithm<Eigen::VectorXd, Eigen::VectorXd, double, T4>::Algorithm(
              algorithm_type, model_type, max_iter, primary_model_fit_max_iter, primary_model_fit_epsilon, warm_start,
              exchange_num, always_select, splicing_type, sub_search){};

    ~abessDAG(){};

    int get_beta_size(int n, int p) { return p * p; }

    Eigen::VectorXi inital_screening(T4 &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, double &coef0,
                                     Eigen::VectorXi &A, Eigen::VectorXi &I, Eigen::VectorXd &bd,
                                     Eigen::VectorXd &weights, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size,
                                     int &N) {
        if (bd.size() == 0) {
            // variable initialization
            int p = X.cols();
            bd = Eigen::VectorXd::Zero(N);

            // pearson correlation
            MatrixXd X1 = X;
            MatrixXd centered = X1.rowwise() - X1.colwise().mean();
            MatrixXd Sigma = (centered.adjoint() * centered);
            VectorXd sd = Sigma.diagonal().eval().cwiseSqrt();
            for (int j = 0; j < p; j++) {
                for (int i = 0; i < j; i++) {
                    bd(j * p + i) = Sigma(i, j) / sd(i) / sd(j);
                }
            }

            // diag
            for (int i = 0; i < p; i++) {
                bd(i * p + i) = -DBL_MAX;
            }

            // alway_select
            for (int i = 0; i < this->always_select.size(); i++) {
                bd(this->always_select(i)) = DBL_MAX;
            }
            // A_init
            for (int i = 0; i < A.size(); i++) {
                bd(A(i)) = DBL_MAX;
            }
        }

        // get Active-set A according to max_k bd
        Eigen::VectorXi A_new = max_k(bd, this->sparsity_level);

        return A_new;
    }

    bool primary_model_fit(T4 &x, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                           double loss0, Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size) {
        cout << "  --> primary fit | A = " << A.transpose() << endl;
        int n = x.rows();
        int p = x.cols();

        cout << "    Test acyclic:" << endl;
        VectorXd temp = VectorXd::Ones(A.size());
        MatrixXd Adj = compute_Adj(temp, A, p);
        if (this->is_cyclic(Adj)) {
            cout << "      Cyclic!" << endl;
            return false;
        }
        cout << "      Pass." << endl;

        // update parents
        vector<VectorXi> parents = this->compute_parents(A, p);

        // compute beta
        int ind = 0;
        for (int i = 0; i < p; i++) {
            int len = parents[i].size();
            if (len != 0) {
                T4 Xpar = X_seg(x, n, parents[i], 0);
                VectorXd y = x.col(i);
                beta.segment(ind, len) = this->lm(Xpar, y);
                ind += len;
            }
        }

        // cout << "    beta_new: " << beta.transpose() << endl;
        return true;
    };

    double loss_function(T4 &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                         Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size, double lambda) {
        MatrixXd Adj = this->compute_Adj(beta, A, X.cols());
        double loss = (X - X * Adj).squaredNorm();
        cout << "  --> loss = " << loss << endl;
        return loss;
    };

    void sacrifice(T4 &X, T4 &XA, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &beta_A, double &coef0,
                   Eigen::VectorXi &A, Eigen::VectorXi &I, Eigen::VectorXd &weights, Eigen::VectorXi &g_index,
                   Eigen::VectorXi &g_size, int N, Eigen::VectorXi &A_ind, Eigen::VectorXd &bd, Eigen::VectorXi &U,
                   Eigen::VectorXi &U_ind, int num) {
        cout << "  --> sacrifice" << endl;
        int n = X.rows();
        int p = X.cols();
        MatrixXd Adj = this->compute_Adj(beta_A, A, p);
        vector<VectorXi> parents = this->compute_parents(A, p);

        // backward
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            // drop A(i) and refit
            VectorXi par(parents[mj].size() - 1);
            int ind = 0;
            for (int k = 0; k < parents[mj].size(); k++) {
                if (parents[mj](k) != mi) par(ind++) = parents[mj](k);
            }
            T4 Xpar = X_seg(X, n, par, 0);
            VectorXd y = X.col(mj);
            VectorXd est_new = Xpar * this->lm(Xpar, y);
            VectorXd est_old = X * Adj.col(mj);
            // loss change
            bd(A(i)) = (est_new - y).squaredNorm() - (est_old - y).squaredNorm();
            // bd(A(i)) = ((est_new + est_old - 2 * y).array() * (est_new - est_old).array()).sum();
            // cout << "    backward: (" << A(i) << ") = " << bd(A(i)) << endl;
        }

        // forward
        for (int i = 0; i < I.size(); i++) {
            int mi = I(i) % p;
            int mj = int(I(i) / p);
            if (mi == mj) {
                bd(I(i)) = -DBL_MAX;
                continue;
            }
            // test acylic
            VectorXd temp = VectorXd::Ones(A.size() + 1);
            VectorXi A_temp(A.size() + 1);
            A_temp.head(A.size()) = A;
            A_temp(A.size()) = I(i);
            MatrixXd Adj_temp = compute_Adj(temp, A_temp, p);
            if (this->is_cyclic(Adj_temp)) {
                bd(I(i)) = -DBL_MAX;
                continue;
            }
            // add I(i) and refit
            VectorXi par(parents[mj].size() + 1);
            par.head(parents[mj].size()) = parents[mj];
            par(parents[mj].size()) = mi;
            T4 Xpar = X_seg(X, n, par, 0);
            VectorXd y = X.col(mj);
            VectorXd est_new = Xpar * this->lm(Xpar, y);
            VectorXd est_old = X * Adj.col(mj);
            // loss change
            bd(I(i)) = (est_old - y).squaredNorm() - (est_new - y).squaredNorm();
            // bd(I(i)) = ((est_old + est_new - 2 * y).array() * (est_old - est_new).array()).sum();
        }

        // // forward
        // MatrixXd E = X - X * Adj;
        // E = E.rowwise() - E.colwise().mean();
        // MatrixXd Omega = (E.colwise().norm()).asDiagonal();
        // MatrixXd inv = (MatrixXd::Identity(p, p) - Adj).inverse();  // TODO(hjh): sparse inverse?
        // MatrixXd temp = Omega - inv.adjoint() * this->Sigma * inv;
        // for (int i = 0; i < A.size(); i++) {
        //     int mi = A(i) % p;
        //     int mj = int(A(i) / p);
        //     bd(A(i)) = abs(temp(mi, mj));
        // }

        // // backward
        // for (int i = 0; i < I.size(); i++) {
        //     int mi = I(i) % p;
        //     int mj = int(I(i) / p);
        //     bd(I(i)) = abs(Adj(mi, mj));
        // }

        // // diag
        // for (int i = 0; i < p; i++) {
        //     bd(i * p + i) = -DBL_MAX;
        // }
    };

    MatrixXd compute_Adj(VectorXd &beta_A, VectorXi &A, int p) {
        MatrixXd Adj = MatrixXd::Zero(p, p);
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            Adj(mi, mj) = beta_A(i);
        }
        return Adj;
    }

    vector<VectorXi> compute_parents(VectorXi &A, int p) {
        VectorXi temp(p);
        vector<VectorXi> parents(p);
        int len = 0, col = 0;
        for (int i = 0; i < A.size(); i++) {
            while (col != int(A(i) / p)) {  // next node
                parents[col++] = temp.segment(0, len).eval();
                len = 0;
            }
            temp(len++) = A(i) % p;
        }

        parents[col++] = temp.segment(0, len).eval();
        while (col < p) {
            parents[col++] = VectorXi::Zero(0);
        }

        return parents;
    }

    VectorXd lm(T4 &XA, VectorXd &y) {
        MatrixXd XTX = XA.adjoint() * XA;
        VectorXd XTY = XA.adjoint() * y;
        return XTX.ldlt().solve(XTY);
    }

    bool is_cyclic(MatrixXd Adj) {
        Adj = Adj.cwiseAbs().eval();
        VectorXd node_in = Adj.colwise().sum();

        int num = node_in.size();
        while (num > 1) {
            // search for "node_in > 0"
            VectorXi ind(num);
            int new_num = 0;
            for (int i = 0; i < num; i++) {
                if (node_in(i) > 0) ind(new_num++) = i;
            }
            if (new_num == num) return true;
            if (new_num == 0) break;
            // delete "node_in == 0"
            Adj = slice_Adj(Adj, ind.head(new_num));
            node_in = Adj.colwise().sum();
            num = node_in.size();
        }
        return false;
    }

    MatrixXd slice_Adj(MatrixXd &Adj, VectorXi ind) {
        int isize = ind.size();
        MatrixXd Adj_new(isize, isize);
        for (int i = 0; i < isize; i++) {
            for (int j = 0; j < isize; j++) {
                Adj_new(i, j) = Adj(ind(i), ind(j));
            }
        }
        return Adj_new;
    }
};

#endif  // SRC_ALGORITHMDAG_H
