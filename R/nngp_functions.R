# Assicurati di eseguire: Rcpp::sourceCpp("../R/nngp_cpp.cpp")

build_nngp_neighbors <- function(coords, m = 15) {
  N <- nrow(coords)

  # Ordinamento spaziale
  ord <- order(coords[, 1], coords[, 2])
  coords_ord <- coords[ord, ]

  cat("Costruzione struttura NNGP (Solo Indici) in C++ per", N, "vertici (m =", m, ")...\n")
  cpp_res <- build_nngp_neighbors_cpp(coords_ord, m)

  return(list(
    coords_ord = coords_ord,
    ord_indices = ord,
    nn_idx = cpp_res$nn_idx,
    m = m
  ))
}

build_nngp_precision <- function(nngp_struct, sigma2, phi) {
  N <- nrow(nngp_struct$coords_ord)

  # Generazione Vettori in C++ (On-the-fly)
  cpp_results <- build_precision_vectors_cpp(
    N = as.integer(N),
    coords_ord = as.matrix(nngp_struct$coords_ord),
    nn_idx = nngp_struct$nn_idx,
    sigma2 = as.numeric(sigma2),
    phi = as.numeric(phi)
  )

  # Assemblaggio Sparso: B_scaled = F^{-1/2} * B_nngp
  B_nngp <- Matrix::sparseMatrix(
    i = cpp_results$i_vec, j = cpp_results$j_vec, x = cpp_results$x_vec, dims = c(N, N)
  )
  B_scaled <- Matrix::Diagonal(x = 1 / sqrt(cpp_results$F_vec)) %*% B_nngp

  # crossprod restituisce direttamente B^T * F^{-1} * B
  Q_nngp <- Matrix::crossprod(B_scaled)

  ord <- nngp_struct$ord_indices
  inv_ord <- order(ord)

  # Forziamo la classe simmetrica dsCMatrix dopo il riordinamento per velocizzare chol()
  Q_out <- Matrix::forceSymmetric(Q_nngp[inv_ord, inv_ord], uplo = "U")

  # Log-determinante analitico O(N) gratuito: log|det(Q)| = -sum(log(F_i))
  attr(Q_out, "log_det") <- -sum(log(cpp_results$F_vec))

  return(Q_out)
}
