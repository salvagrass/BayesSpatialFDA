#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// MODULO 1: Ricerca Vicini e Caching delle Distanze (Sostituisce build_nngp_neighbors)
// [[Rcpp::export]]
List build_nngp_neighbors_cpp(arma::mat coords_ord, int m) {
  int N = coords_ord.n_rows;

  // Matrici per indici e distanze con NA
  IntegerMatrix nn_idx(N, m);
  NumericMatrix nn_dist(N, m);
  std::fill(nn_idx.begin(), nn_idx.end(), NA_INTEGER);
  std::fill(nn_dist.begin(), nn_dist.end(), NA_REAL);

  List nn_dist_matrix(N);
  nn_dist_matrix[0] = arma::mat(1, 1, arma::fill::zeros);

  for (int i = 1; i < N; i++) {
    int num_neighbors = std::min(m, i);
    arma::vec dists(i);

    // Calcolo distanze dal punto i a tutti i punti precedenti
    for (int j = 0; j < i; j++) {
      dists[j] = std::sqrt(std::pow(coords_ord(j, 0) - coords_ord(i, 0), 2) +
        std::pow(coords_ord(j, 1) - coords_ord(i, 1), 2));
    }

    // Trova i vicini più prossimi
    arma::uvec sort_idx = arma::sort_index(dists);
    arma::uvec nearest = sort_idx.head(num_neighbors);

    arma::mat coords_neigh(num_neighbors, 2);
    for (int k = 0; k < num_neighbors; k++) {
      int n_idx = nearest[k];
      nn_idx(i, k) = n_idx + 1; // +1 per indici 1-based di R
      nn_dist(i, k) = dists[n_idx];
      coords_neigh(k, 0) = coords_ord(n_idx, 0);
      coords_neigh(k, 1) = coords_ord(n_idx, 1);
    }

    // Caching delle distanze reciproche TRA i vicini
    arma::mat dist_mat(num_neighbors, num_neighbors, arma::fill::zeros);
    if (num_neighbors > 1) {
      for (int r = 0; r < num_neighbors; r++) {
        for (int c = r + 1; c < num_neighbors; c++) {
          double d = std::sqrt(std::pow(coords_neigh(r, 0) - coords_neigh(c, 0), 2) +
                               std::pow(coords_neigh(r, 1) - coords_neigh(c, 1), 2));
          dist_mat(r, c) = d;
          dist_mat(c, r) = d;
        }
      }
    }
    nn_dist_matrix[i] = dist_mat;
  }

  return List::create(
    Named("nn_idx") = nn_idx,
    Named("nn_dist") = nn_dist,
    Named("nn_dist_matrix") = nn_dist_matrix
  );
}

// MODULO 2: Costruzione Veloce della Matrice di Precisione
// [[Rcpp::export]]
List build_precision_vectors_cpp(int N, IntegerMatrix nn_idx, List nn_dist_matrix,
                                 NumericMatrix nn_dist, double sigma2, double phi) {

  std::vector<int> i_vec; std::vector<int> j_vec; std::vector<double> x_vec;
  i_vec.reserve(N * 15); j_vec.reserve(N * 15); x_vec.reserve(N * 15);

  NumericVector F_vec(N);
  F_vec[0] = sigma2;
  i_vec.push_back(1); j_vec.push_back(1); x_vec.push_back(1.0);

  for (int i = 1; i < N; i++) {
    int num_neigh = 0;
    for (int k = 0; k < nn_idx.ncol(); k++) {
      if (!IntegerVector::is_na(nn_idx(i, k))) num_neigh++;
      else break;
    }

    arma::uvec neigh(num_neigh);
    arma::vec dist_i_neigh(num_neigh);
    for (int k = 0; k < num_neigh; k++) {
      neigh[k] = nn_idx(i, k);
      dist_i_neigh[k] = nn_dist(i, k);
    }

    NumericMatrix dist_mat_R = nn_dist_matrix[i];
    arma::mat dist_neigh = as<arma::mat>(dist_mat_R);

    arma::mat C_neigh = sigma2 * arma::exp(-phi * dist_neigh);
    arma::vec c_cross = sigma2 * arma::exp(-phi * dist_i_neigh);

    arma::vec B_i = arma::solve(C_neigh, c_cross);
    F_vec[i] = sigma2 - arma::dot(c_cross, B_i);

    i_vec.push_back(i + 1); j_vec.push_back(i + 1); x_vec.push_back(1.0);
    for(int k = 0; k < num_neigh; k++) {
      i_vec.push_back(i + 1); j_vec.push_back(neigh[k]); x_vec.push_back(-B_i[k]);
    }
  }

  return List::create(Named("i_vec") = i_vec, Named("j_vec") = j_vec,
                      Named("x_vec") = x_vec, Named("F_vec") = F_vec);
}
