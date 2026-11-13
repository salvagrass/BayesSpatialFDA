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
    ord_indices = as.integer(ord),
    nn_idx = cpp_res$nn_idx,
    m = m
  ))
}

# Assembla la matrice sparsa solo quando necessario (inizializzazione o accettazione Metropolis)
assemble_nngp_precision <- function(cpp_results, N) {
  B_scaled <- Matrix::sparseMatrix(
    i = cpp_results$i_vec,
    j = cpp_results$j_vec,
    x = cpp_results$x_vec,
    dims = c(N, N)
  )
  # Restituisce dgCMatrix già nell'ordine originale dei nodi (zero permutazioni inv_ord)
  Q_out <- Matrix::t(B_scaled) %*% B_scaled
  attr(Q_out, "log_det") <- cpp_results$log_det
  return(Q_out)
}

build_nngp_precision <- function(nngp_struct, sigma2, phi) {
  N <- nrow(nngp_struct$coords_ord)

  cpp_results <- build_precision_vectors_cpp(
    N = as.integer(N),
    coords_ord = as.matrix(nngp_struct$coords_ord),
    nn_idx = nngp_struct$nn_idx,
    ord = nngp_struct$ord_indices,
    sigma2 = as.numeric(sigma2),
    phi = as.numeric(phi)
  )

  return(assemble_nngp_precision(cpp_results, N))
}
