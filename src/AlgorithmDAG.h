#ifndef SRC_ALGORITHMDAG_H
#define SRC_ALGORITHMDAG_H

#include <Spectra/SymEigsSolver.h>

#include "Algorithm.h"

using namespace Spectra;

template <class T4>
class abessDAG : public Algorithm<Eigen::VectorXd, Eigen::VectorXd, double, T4> {
   public:
    // MatrixXd Adj;  // adjacency matrix

    abessDAG(int algorithm_type, int model_type, int max_iter = 30, int primary_model_fit_max_iter = 10,
             double primary_model_fit_epsilon = 1e-8, bool warm_start = true, int exchange_num = 5,
             Eigen::VectorXi always_select = Eigen::VectorXi::Zero(0), int splicing_type = 1, int sub_search = 0)
        : Algorithm<Eigen::VectorXd, Eigen::VectorXd, double, T4>::Algorithm(
              algorithm_type, model_type, max_iter, primary_model_fit_max_iter, primary_model_fit_epsilon, warm_start,
              exchange_num, always_select, splicing_type, sub_search){};

    ~abessDAG(){};

    // void inital_setting(T4 &X, VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXi &g_index, Eigen::VectorXi
    // &g_size,
    //                     int &N) {

    // }

    bool primary_model_fit(T4 &x, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                           double loss0, Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size) {
        int n = x.rows();
        int p = x.cols();
        MatrixXd Adj = this->compute_Adj(beta, A, n, p);

        beta = this->compute_beta(Adj, p, A.size());
    };

    double loss_function(T4 &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0,
                         Eigen::VectorXi &A, Eigen::VectorXi &g_index, Eigen::VectorXi &g_size, double lambda) {
        int n = x.rows();
        int p = x.cols();
        MatrixXd Adj = this->compute_Adj(beta, A, n, p);
        MatrixXd Xj(n, 1);
        MatrixXd Xother(n, p - 1);
        for (int j = 0; j < p; j++){
            X_est = X * Adj;

        }
    };

    void sacrifice(T4 &X, T4 &XA, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &beta_A, double &coef0,
                   Eigen::VectorXi &A, Eigen::VectorXi &I, Eigen::VectorXd &weights, Eigen::VectorXi &g_index,
                   Eigen::VectorXi &g_size, int N, Eigen::VectorXi &A_ind, Eigen::VectorXd &bd, Eigen::VectorXi &U,
                   Eigen::VectorXi &U_ind, int num) {
        int n = x.rows();
        int p = x.cols();
        MatrixXd Adj = this->compute_Adj(beta_A, A, n, p);
    };

    MatrixXd compute_Adj(VectorXd &beta_A, VectorXi &A, int p) {
        MatrixXd Adj = MatrixXd::Zero(p, p);
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            Adj(mi, mj) = beta(i);
        }
        return Adj;
    }

    VectorXd compute_beta_A(MatrixXd &Adj, VectorXi &A, int beta_A_size) {
        VectorXd beta_A = VectorXd::Zero(beta_A_size);
        int p = Adj.rows();
        for (int i = 0; i < A.size(); i++) {
            int mi = A(i) % p;
            int mj = int(A(i) / p);
            beta_A(i) = Adj(mi, mj);
        }
        return beta_A;
    }

    void split_col(MatrixXd &X, int j, MatrixXd &Xj, MatrixXd &Xother){
        Xj.col(0) = X.col(j);
        Xother << X.block(0, 0, X.rows(), j), X.block(0, j + 1, X.rows(), X.cols()-j-1);
    }
};

#endif  // SRC_ALGORITHMDAG_H
