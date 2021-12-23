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

    void inital_setting(T4 &X, VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size,
                        int &N) {
        MatrixXd X1 = X;
        MatrixXd centered = X1.rowwise() - X1.colwise().mean();
        this->Sigma = (centered.adjoint() * centered);
    }

    bool primary_model_fit(T4 &x, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                           double loss0, Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size) {
        int n = x.rows();
        int p = x.cols();
        // MatrixXd Adj = this->compute_Adj(beta, A, n, p);

        MatrixXd X_temp(n, p);
        int ind = 0, col = -1, beta_ind = 0;
        for (int i = 0; i < A.size(); i++) {
            if (col != int(A(i) / p)) {  // new node
                if (ind > 0) {           // update beta
                    beta.segment(beta_ind, ind) = lm(X_temp.block(0, 0, n, ind), x.col(col));
                    beta_ind += ind;
                }
                col = int(A(i) / p);
                ind = 0;
            }
            X_temp.col(ind++) = x.col(A(i) % p);
        }
    };

    double loss_function(T4 &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                         Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size, double lambda) {
        int n = X.rows();
        int p = X.cols();
        MatrixXd Adj = this->compute_Adj(beta, A, p);
        return (X - X * Adj).norm();
    };

    void sacrifice(T4 &X, T4 &XA, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &beta_A, double &coef0,
                   Eigen::VectorXi &A, Eigen::VectorXi &I, Eigen::VectorXd &weights, Eigen::VectorXi &g_index,
                   Eigen::VectorXi &g_size, int N, Eigen::VectorXi &A_ind, Eigen::VectorXd &bd, Eigen::VectorXi &U,
                   Eigen::VectorXi &U_ind, int num) {
        int n = X.rows();
        int p = X.cols();
        MatrixXd Adj = this->compute_Adj(beta_A, A, p);

        // forward
        MatrixXd E = X - X * Adj;
        E = E.rowwise() - E.colwise().mean();
        MatrixXd Omega = (E.colwise().squaredNorm()).asDiagonal();
        MatrixXd inv = (MatrixXd::Identity(p, p) - Adj).inverse();  // TODO(hjh): sparse inverse?
        MatrixXd temp = Omega - inv.adjoint() * this->Sigma * inv;
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            bd(A(i)) = abs(temp(mi, mj));
        }

        // backward
        for (int i = 0; i < I.size(); i++) {
            int mi = I(i) % p;
            int mj = int(I(i) / p);
            bd(I(i)) = abs(Adj(mi, mj));
        }
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

    // VectorXd compute_beta_A(MatrixXd &Adj, VectorXi &A, int beta_A_size) {
    //     VectorXd beta_A = VectorXd::Zero(beta_A_size);
    //     int p = Adj.rows();
    //     for (int i = 0; i < A.size(); i++) {
    //         int mi = A(i) % p;
    //         int mj = int(A(i) / p);
    //         beta_A(i) = Adj(mi, mj);
    //     }
    //     return beta_A;
    // }

    VectorXd lm(MatrixXd X, VectorXd y) { return (X.adjoint() * X).ldlt().solve(X.adjoint() * y); }
};

#endif  // SRC_ALGORITHMDAG_H
