# Assicurati che RcppArmadillo sia installato e il file sia compilato:
# Rcpp::sourceCpp("nngp_cpp.cpp")

build_nngp_neighbors <- function(coords, m = 15) {
  N <- nrow(coords)

  # 1. Ordinamento spaziale (rimane in R per semplicità)
  ord <- order(coords[, 1], coords[, 2])
  coords_ord <- coords[ord, ]

  # 2. Computazione massiva in C++
  cat("Costruzione struttura NNGP in C++ per", N, "vertici (m =", m, ")...\n")
  cpp_res <- build_nngp_neighbors_cpp(coords_ord, m)

  return(list(
    coords_ord = coords_ord,
    ord_indices = ord,
    nn_idx = cpp_res$nn_idx,
    nn_dist = cpp_res$nn_dist,
    nn_dist_matrix = cpp_res$nn_dist_matrix,
    m = m
  ))
}

build_nngp_precision <- function(nngp_struct, sigma2, phi) {
  N <- nrow(nngp_struct$coords_ord)

  # 1. Generazione Vettori in C++
  cpp_results <- build_precision_vectors_cpp(
    N = as.integer(N),
    nn_idx = nngp_struct$nn_idx,
    nn_dist_matrix = nngp_struct$nn_dist_matrix,
    nn_dist = nngp_struct$nn_dist,
    sigma2 = as.numeric(sigma2),
    phi = as.numeric(phi)
  )

  # 2. Assemblaggio Sparso in R (Sfrutta C-SuiteSparse)
  B_nngp <- Matrix::sparseMatrix(
    i = cpp_results$i_vec, j = cpp_results$j_vec, x = cpp_results$x_vec, dims = c(N, N)
  )
  F_inv <- Matrix::Diagonal(x = 1 / cpp_results$F_vec)
  Q_nngp <- Matrix::t(B_nngp) %*% F_inv %*% B_nngp

  ord <- nngp_struct$ord_indices
  inv_ord <- order(ord)
  return(Q_nngp[inv_ord, inv_ord])
}
