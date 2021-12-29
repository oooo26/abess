#ifndef SRC_ALGORITHMDAG_H
#define SRC_ALGORITHMDAG_H

#include <Spectra/SymEigsSolver.h>

#include "Algorithm.h"

using namespace Spectra;

template <class T4>
class abessDAG : public Algorithm<Eigen::VectorXd, Eigen::VectorXd, double, T4> {
   public:
    // MatrixXd Adj;  // adjacency matrix
    MatrixXd Sigma;

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

            // initial acyclic
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < i; j++) {
                    bd(i * p + j) = 1;
                }
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

    void inital_setting(T4 &X, VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size,
                        int &N) {
        MatrixXd X1 = X;
        MatrixXd centered = X1.rowwise() - X1.colwise().mean();
        this->Sigma = (centered.adjoint() * centered);
    }

    bool primary_model_fit(T4 &x, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                           double loss0, Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size) {
        cout << "  --> primary fit | A = " << A.transpose() << endl;
        int n = x.rows();
        int p = x.cols();
        // MatrixXd Adj = this->compute_Adj(beta, A, p);

        if (this->is_cyclic(A, p)) {
            cout << "    cyclic !" << endl;
            return false;
        }

        MatrixXd X_temp(n, p);
        int ind = 0, col = -1, beta_ind = 0;
        for (int i = 0; i < A.size(); i++) {
            if (col != int(A(i) / p)) {  // new node
                if (ind > 0) {           // update beta
                    beta.segment(beta_ind, ind) = this->lm(X_temp.block(0, 0, n, ind), x.col(col));
                    beta_ind += ind;
                }
                col = int(A(i) / p);
                ind = 0;
            }
            X_temp.col(ind++) = x.col(A(i) % p);
        }
        beta.segment(beta_ind, ind) = this->lm(X_temp.block(0, 0, n, ind), x.col(col));
        return true;
    };

    double loss_function(T4 &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                         Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size, double lambda) {
        cout << "  --> loss" << endl;
        // int n = X.rows();
        int p = X.cols();
        MatrixXd Adj = this->compute_Adj(beta, A, p);
        cout << "    " << (X - X * Adj).norm() << endl;
        return (X - X * Adj).norm();
    };

    void sacrifice(T4 &X, T4 &XA, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &beta_A, double &coef0,
                   Eigen::VectorXi &A, Eigen::VectorXi &I, Eigen::VectorXd &weights, Eigen::VectorXi &g_index,
                   Eigen::VectorXi &g_size, int N, Eigen::VectorXi &A_ind, Eigen::VectorXd &bd, Eigen::VectorXi &U,
                   Eigen::VectorXi &U_ind, int num) {
        cout << "  --> sacrifice" << endl;
        int n = X.rows();
        int p = X.cols();
        MatrixXd Adj = this->compute_Adj(beta_A, A, p);

        // backward
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            VectorXd est_old = X * Adj.col(mj).eval();
            VectorXd est_new = est_old - beta_A(i) * X.col(mi).eval();
            VectorXd Xj = X.col(mj);
            bd(A(i)) = (est_new - Xj).squaredNorm() - (est_old - Xj).squaredNorm();
            cout<<"    backward: ("<<A(i)<<") = "<<bd(A(i))<<endl;
            // bd(A(i)) = ((est_new + est_old - 2 * Xj).array() * (est_new - est_old).array()).sum();
        }

        // forward
        for (int i = 0; i < I.size(); i++) {
            int mi = I(i) % p;
            int mj = int(I(i) / p);
            if (mi == mj) {
                bd(I(i)) = -DBL_MAX;
                continue;
            }
            VectorXd est_old = X * Adj.col(mj).eval();
            VectorXd est_new(p);
            VectorXd Xj = X.col(mj);
            MatrixXd X_temp(n, p);
            int ind = 1;
            X_temp.col(0) = X.col(mi);
            for (int k = 0; k < p; k++) {
                if (Adj(k, mj) != 0) {
                    X_temp.col(ind++) = X.col(k);
                }
            }
            VectorXd beta_temp = this->lm(X_temp.block(0, 0, n, ind), X.col(mj));
            ind = 1;
            est_new(mi) = beta_temp(0);
            for (int k = 0; k < p; k++) {
                if (Adj(k, mj) != 0) {
                    est_new(k) = beta_temp(ind++);
                }
            }
            bd(I(i)) = (est_new - Xj).squaredNorm() - (est_old - Xj).squaredNorm();
            // bd(I(i)) = ((est_new + est_old - 2 * Xj).array() * (est_new - est_old).array()).sum();
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

    VectorXd lm(MatrixXd X, VectorXd y) { return (X.adjoint() * X).ldlt().solve(X.adjoint() * y); }

    bool is_cyclic(VectorXi &A, int p) {
        VectorXd temp = VectorXd::Ones(A.size());
        MatrixXd Adj = this->compute_Adj(temp, A, p);
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
