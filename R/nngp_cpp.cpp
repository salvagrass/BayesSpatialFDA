#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Cache statica per la fattorizzazione simbolica del grafo della mesh (Step 2)
static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> g_llt_solver;
static bool g_pattern_initialized = false;
static int g_cached_N = -1;
static int g_cached_nnz = -1;

// MODULO 1: Ricerca Vicini On-the-fly (Ottimizzata con partial_sort)
// [[Rcpp::export]]
List build_nngp_neighbors_cpp(MatrixXd coords_ord, int m) {
  // Resetta la cache simbolica quando viene costruita una nuova mesh
  g_pattern_initialized = false;

  int N = coords_ord.rows();
  IntegerMatrix nn_idx(N, m);
  std::fill(nn_idx.begin(), nn_idx.end(), NA_INTEGER);

  for (int i = 1; i < N; i++) {
    int num_neighbors = std::min(m, i);
    std::vector<double> dists_sq(i);

    for (int j = 0; j < i; j++) {
      double dx = coords_ord(j, 0) - coords_ord(i, 0);
      double dy = coords_ord(j, 1) - coords_ord(i, 1);
      dists_sq[j] = dx * dx + dy * dy;
    }

    std::vector<int> sort_idx(i);
    std::iota(sort_idx.begin(), sort_idx.end(), 0);
    std::partial_sort(sort_idx.begin(), sort_idx.begin() + num_neighbors, sort_idx.end(),
                      [&dists_sq](int a, int b) {
                        return dists_sq[a] < dists_sq[b];
                      });

    for (int k = 0; k < num_neighbors; k++) {
      nn_idx(i, k) = sort_idx[k] + 1;
    }
  }
  return List::create(Named("nn_idx") = nn_idx);
}

// MODULO 2: Costruzione Vettori B_scaled = F^{-1/2} * B già mappati nell'ordine originale
// [[Rcpp::export]]
List build_precision_vectors_cpp(int N, MatrixXd coords_ord, IntegerMatrix nn_idx,
                                 IntegerVector ord, double sigma2, double phi) {

  int max_nn = nn_idx.ncol();
  std::vector<int> nn_counts(N, 0);
  std::vector<int> start_idx(N, 0);

  std::vector<int> nn_idx_std(N * max_nn, -1);
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < max_nn; k++) {
      if (!IntegerVector::is_na(nn_idx(i, k))) {
        nn_idx_std[i * max_nn + k] = nn_idx(i, k);
      }
    }
  }

  int total_elements = 1;
  for (int i = 1; i < N; i++) {
    int c = 0;
    for (int k = 0; k < max_nn; k++) {
      if (nn_idx_std[i * max_nn + k] != -1) c++;
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
  double inv_sqrt_F0 = 1.0 / std::sqrt(F_vec[0]);
  i_vec[0] = 1;
  j_vec[0] = ord[0];
  x_vec[0] = inv_sqrt_F0;
  double log_det = -std::log(F_vec[0]);

  for (int i = 1; i < N; i++) {
    int c = nn_counts[i];
    int s = start_idx[i];

    MatrixXd coords_neigh(c, 2);
    VectorXd dist_i_neigh(c);

    for (int k = 0; k < c; k++) {
      int n_idx = nn_idx_std[i * max_nn + k] - 1;
      coords_neigh(k, 0) = coords_ord(n_idx, 0);
      coords_neigh(k, 1) = coords_ord(n_idx, 1);

      double dx = coords_ord(i, 0) - coords_ord(n_idx, 0);
      double dy = coords_ord(i, 1) - coords_ord(n_idx, 1);
      dist_i_neigh[k] = std::sqrt(dx * dx + dy * dy);
    }

    MatrixXd dist_neigh = MatrixXd::Zero(c, c);
    for (int r = 0; r < c; r++) {
      for (int col = r + 1; col < c; col++) {
        double dx = coords_neigh(r, 0) - coords_neigh(col, 0);
        double dy = coords_neigh(r, 1) - coords_neigh(col, 1);
        double d = std::sqrt(dx * dx + dy * dy);
        dist_neigh(r, col) = d;
        dist_neigh(col, r) = d;
      }
    }

    MatrixXd C_neigh = sigma2 * (-phi * dist_neigh.array()).exp();
    VectorXd c_cross = sigma2 * (-phi * dist_i_neigh.array()).exp();

    VectorXd B_i = C_neigh.llt().solve(c_cross);
    F_vec[i] = sigma2 - c_cross.dot(B_i);
    double inv_sqrt_Fi = 1.0 / std::sqrt(F_vec[i]);
    log_det -= std::log(F_vec[i]);

    i_vec[s] = i + 1;
    j_vec[s] = ord[i];
    x_vec[s] = inv_sqrt_Fi;

    for (int k = 0; k < c; k++) {
      int n_idx = nn_idx_std[i * max_nn + k] - 1;
      i_vec[s + 1 + k] = i + 1;
      j_vec[s + 1 + k] = ord[n_idx];
      x_vec[s + 1 + k] = -B_i[k] * inv_sqrt_Fi;
    }
  }

  return List::create(Named("i_vec") = i_vec, Named("j_vec") = j_vec,
                      Named("x_vec") = x_vec, Named("F_vec") = F_vec,
                      Named("log_det") = log_det);
}

// MODULO 3: Valutazione Lazy di phi_prop (Step 5) senza costruire matrici sparse in R
// [[Rcpp::export]]
List eval_nngp_proposal_cpp(int N, MatrixXd coords_ord, IntegerMatrix nn_idx,
                            IntegerVector ord, double sigma2, double phi,
                            const MatrixXd& beta) {

  int max_nn = nn_idx.ncol();
  std::vector<int> nn_counts(N, 0);
  std::vector<int> start_idx(N, 0);

  std::vector<int> nn_idx_std(N * max_nn, -1);
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < max_nn; k++) {
      if (!IntegerVector::is_na(nn_idx(i, k))) {
        nn_idx_std[i * max_nn + k] = nn_idx(i, k);
      }
    }
  }

  int total_elements = 1;
  for (int i = 1; i < N; i++) {
    int c = 0;
    for (int k = 0; k < max_nn; k++) {
      if (nn_idx_std[i * max_nn + k] != -1) c++;
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
  double inv_sqrt_F0 = 1.0 / std::sqrt(F_vec[0]);
  int orig_0 = ord[0];
  i_vec[0] = 1;
  j_vec[0] = orig_0;
  x_vec[0] = inv_sqrt_F0;

  double log_det = -std::log(F_vec[0]);
  double quad_form = beta.row(orig_0 - 1).squaredNorm() / F_vec[0];

  for (int i = 1; i < N; i++) {
    int c = nn_counts[i];
    int s = start_idx[i];
    int orig_i = ord[i];

    MatrixXd coords_neigh(c, 2);
    VectorXd dist_i_neigh(c);

    for (int k = 0; k < c; k++) {
      int n_idx = nn_idx_std[i * max_nn + k] - 1;
      coords_neigh(k, 0) = coords_ord(n_idx, 0);
      coords_neigh(k, 1) = coords_ord(n_idx, 1);

      double dx = coords_ord(i, 0) - coords_ord(n_idx, 0);
      double dy = coords_ord(i, 1) - coords_ord(n_idx, 1);
      dist_i_neigh[k] = std::sqrt(dx * dx + dy * dy);
    }

    MatrixXd dist_neigh = MatrixXd::Zero(c, c);
    for (int r = 0; r < c; r++) {
      for (int col = r + 1; col < c; col++) {
        double dx = coords_neigh(r, 0) - coords_neigh(col, 0);
        double dy = coords_neigh(r, 1) - coords_neigh(col, 1);
        double d = std::sqrt(dx * dx + dy * dy);
        dist_neigh(r, col) = d;
        dist_neigh(col, r) = d;
      }
    }

    MatrixXd C_neigh = sigma2 * (-phi * dist_neigh.array()).exp();
    VectorXd c_cross = sigma2 * (-phi * dist_i_neigh.array()).exp();

    VectorXd B_i = C_neigh.llt().solve(c_cross);
    F_vec[i] = sigma2 - c_cross.dot(B_i);
    double inv_sqrt_Fi = 1.0 / std::sqrt(F_vec[i]);

    log_det -= std::log(F_vec[i]);

    RowVectorXd res_i = beta.row(orig_i - 1);
    i_vec[s] = i + 1;
    j_vec[s] = orig_i;
    x_vec[s] = inv_sqrt_Fi;

    for (int k = 0; k < c; k++) {
      int n_idx = nn_idx_std[i * max_nn + k] - 1;
      int orig_n = ord[n_idx];
      i_vec[s + 1 + k] = i + 1;
      j_vec[s + 1 + k] = orig_n;
      x_vec[s + 1 + k] = -B_i[k] * inv_sqrt_Fi;
      res_i.noalias() -= B_i[k] * beta.row(orig_n - 1);
    }

    quad_form += res_i.squaredNorm() / F_vec[i];
  }

  return List::create(Named("i_vec") = i_vec, Named("j_vec") = j_vec,
                      Named("x_vec") = x_vec, Named("F_vec") = F_vec,
                      Named("log_det") = log_det, Named("quad_form") = quad_form);
}

// MODULO 4: Campionamento Veloce Step 2 con Fattorizzazione Simbolica Fissa
// [[Rcpp::export]]
MatrixXd sample_beta_cpp(const Eigen::Map<Eigen::SparseMatrix<double>>& R_phi,
                         double sigma2_spaz,
                         double sigma2_err,
                         const MatrixXd& Y_minus_mu,
                         const MatrixXd& Z_norm) {
  int N = R_phi.rows();
  Eigen::SparseMatrix<double> Q_post = R_phi / sigma2_spaz;

  double inv_err = 1.0 / sigma2_err;
  for (int k = 0; k < N; ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Q_post, k); it; ++it) {
      if (it.row() == k) {
        it.valueRef() += inv_err;
        break;
      }
    }
  }

  // Analisi simbolica del grafo eseguita SOLO alla prima iterazione
  if (!g_pattern_initialized || g_cached_N != N || g_cached_nnz != Q_post.nonZeros()) {
    g_llt_solver.analyzePattern(Q_post);
    g_pattern_initialized = true;
    g_cached_N = N;
    g_cached_nnz = Q_post.nonZeros();
  }

  // Fattorizzazione numerica veloce
  g_llt_solver.factorize(Q_post);

  MatrixXd dest = Y_minus_mu * inv_err;
  if (g_llt_solver.permutationP().size() > 0) {
    dest = g_llt_solver.permutationP() * dest;
  }
  g_llt_solver.matrixL().solveInPlace(dest);
  dest += Z_norm;
  g_llt_solver.matrixU().solveInPlace(dest);
  if (g_llt_solver.permutationPinv().size() > 0) {
    dest = g_llt_solver.permutationPinv() * dest;
  }

  return dest;
}

// MODULO 5: Forma Quadratica Veloce per Step 4
// [[Rcpp::export]]
double quad_form_cpp(const Eigen::Map<Eigen::SparseMatrix<double>>& R_phi,
                     const MatrixXd& beta) {
  return (beta.array() * (R_phi * beta).array()).sum();
}
