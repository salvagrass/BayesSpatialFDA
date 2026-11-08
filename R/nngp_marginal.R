# ==============================================================================
# File: R/nngp_marginal.R
# Funzione per il calcolo della Log-Verosimiglianza Marginale (Collapsed NNGP)
# ==============================================================================

#' Calcola p(Y | mu, sigma2_spaz, phi, sigma2_err) integrando via il campo latente
#'
#' @param sig2_spaz Varianza spaziale proposta
#' @param phi Decadimento spaziale proposto
#' @param sig2_err Varianza dell'errore (rumore bianco)
#' @param Y_mat Matrice dei dati (N x M)
#' @param curr_mu Vettore del campo medio spaziale attuale
#' @param nngp_struct Struttura dei vicini (m=10) pre-calcolata
#' @return Valore scalare della log-verosimiglianza marginale
log_marg_likelihood <- function(sig2_spaz, phi, sig2_err, Y_mat, curr_mu, nngp_struct) {

  N <- nrow(Y_mat)
  M <- ncol(Y_mat)

  # 1. Costruiamo la matrice di precisione sparsa NNGP (Q)
  Q_spaz <- build_nngp_precision(nngp_struct, sig2_spaz, phi)

  # 2. Matrice M per l'Identità di Woodbury
  M_mat <- Q_spaz + Diagonal(N, 1 / sig2_err)
  chol_M <- chol(M_mat) # Cholesky sparsa iper-veloce

  # 3. Calcolo del Determinante Marginale
  ldet_Q <- as.numeric(determinant(Q_spaz, logarithm = TRUE)$modulus)
  ldet_M <- as.numeric(determinant(M_mat, logarithm = TRUE)$modulus)

  # det(Sigma_Y) espresso tramite i determinanti di Q e M
  ldet_SigmaY <- -ldet_Q + N * log(sig2_err) + ldet_M

  # 4. Calcolo della Forma Quadratica Marginale su tutti gli M anni
  quad_sum <- 0
  for(i in 1:M) {
    z_i <- Y_mat[, i] - curr_mu

    # Risolviamo il sistema v = M^(-1) * z in modo efficiente
    v_i <- as.vector(solve(chol_M, solve(t(chol_M), z_i)))

    # Identità di Woodbury per la forma quadratica z^T * Sigma_Y^(-1) * z
    quad_i <- (sum(z_i^2) / sig2_err) - (sum(z_i * v_i) / (sig2_err^2))
    quad_sum <- quad_sum + quad_i
  }

  # Log-Verosimiglianza Gaussiana finale
  return(-0.5 * M * ldet_SigmaY - 0.5 * quad_sum)
}
