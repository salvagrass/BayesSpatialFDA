#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// MODULO 1: Ricerca Vicini On-the-fly
// [[Rcpp::export]]
List build_nngp_neighbors_cpp(MatrixXd coords_ord, int m) {
  int N = coords_ord.rows();
  IntegerMatrix nn_idx(N, m);
  std::fill(nn_idx.begin(), nn_idx.end(), NA_INTEGER);

  for (int i = 1; i < N; i++) {
    int num_neighbors = std::min(m, i);
    std::vector<double> dists(i);

    for (int j = 0; j < i; j++) {
      dists[j] = std::sqrt(std::pow(coords_ord(j, 0) - coords_ord(i, 0), 2) +
        std::pow(coords_ord(j, 1) - coords_ord(i, 1), 2));
    }

    std::vector<int> sort_idx(i);
    std::iota(sort_idx.begin(), sort_idx.end(), 0);
    std::sort(sort_idx.begin(), sort_idx.end(), [&dists](int a, int b) {
      return dists[a] < dists[b];
    });

    for (int k = 0; k < num_neighbors; k++) {
      nn_idx(i, k) = sort_idx[k] + 1;
    }
  }
  return List::create(Named("nn_idx") = nn_idx);
}

// MODULO 2: Costruzione Precisione (Eigen LLT Ultra-Veloce)
// [[Rcpp::export]]
List build_precision_vectors_cpp(int N, MatrixXd coords_ord, IntegerMatrix nn_idx,
                                 double sigma2, double phi) {

  int max_nn = nn_idx.ncol();
  std::vector<int> nn_counts(N, 0);
  std::vector<int> start_idx(N, 0);

  // Trasferiamo la IntegerMatrix di R in un std::vector C++ nativo
  std::vector<int> nn_idx_std(N * max_nn, -1);
  for(int i = 0; i < N; i++) {
    for(int k = 0; k < max_nn; k++) {
      if (!IntegerVector::is_na(nn_idx(i, k))) {
        nn_idx_std[i * max_nn + k] = nn_idx(i, k);
      }
    }
  }

  // Pre-allocazione precisa
  int total_elements = 1;
  for(int i = 1; i < N; i++) {
    int c = 0;
    for(int k = 0; k < max_nn; k++) {
      if(nn_idx_std[i * max_nn + k] != -1) c++;
      else break;
    }
    nn_counts[i] = c;
    start_idx[i] = total_elements;
    total_elements += (1 + c);
  }

  std::vector<int> i_vec(total_elements);
  std::vector<int> j_vec(total_elements);
  std::vector<double> x_vec(total_elements);
  VectorXd F_vec(N);

  F_vec[0] = sigma2;
  i_vec[0] = 1; j_vec[0] = 1; x_vec[0] = 1.0;

  // Ciclo Veloce e Stabile
  for (int i = 1; i < N; i++) {
    int c = nn_counts[i];
    int s = start_idx[i];

    MatrixXd coords_neigh(c, 2);
    VectorXd dist_i_neigh(c);

    for(int k = 0; k < c; k++) {
      int n_idx = nn_idx_std[i * max_nn + k] - 1;
      coords_neigh(k, 0) = coords_ord(n_idx, 0);
      coords_neigh(k, 1) = coords_ord(n_idx, 1);

      dist_i_neigh[k] = std::sqrt(std::pow(coords_ord(i, 0) - coords_ord(n_idx, 0), 2) +
        std::pow(coords_ord(i, 1) - coords_ord(n_idx, 1), 2));
    }

    MatrixXd dist_neigh = MatrixXd::Zero(c, c);
    for(int r = 0; r < c; r++){
      for(int col = r + 1; col < c; col++){
        double d = std::sqrt(std::pow(coords_neigh(r, 0) - coords_neigh(col, 0), 2) +
                             std::pow(coords_neigh(r, 1) - coords_neigh(col, 1), 2));
        dist_neigh(r, col) = d;
        dist_neigh(col, r) = d;
      }
    }

    MatrixXd C_neigh = sigma2 * (-phi * dist_neigh.array()).exp();
    VectorXd c_cross = sigma2 * (-phi * dist_i_neigh.array()).exp();

    VectorXd B_i = C_neigh.llt().solve(c_cross);
    F_vec[i] = sigma2 - c_cross.dot(B_i);

    i_vec[s] = i + 1;
    j_vec[s] = i + 1;
    x_vec[s] = 1.0;

    for(int k = 0; k < c; k++) {
      i_vec[s + 1 + k] = i + 1;
      j_vec[s + 1 + k] = nn_idx_std[i * max_nn + k];
      x_vec[s + 1 + k] = -B_i[k];
    }
  }

  return List::create(Named("i_vec") = i_vec, Named("j_vec") = j_vec,
                      Named("x_vec") = x_vec, Named("F_vec") = F_vec);
}
