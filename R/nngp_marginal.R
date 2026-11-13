# ==============================================================================
# File: R/nngp_marginal.R
# Calcolo Stabile della Verosimiglianza Marginale (Identità di Rue & Held)
# ==============================================================================

log_marg_likelihood_robust <- function(sig2_spaz, phi, sig2_err, Y_mat, curr_mu, nngp_struct) {
  N <- nrow(Y_mat)
  M <- ncol(Y_mat)

  Q_spaz <- build_nngp_precision(nngp_struct, sig2_spaz, phi)
  Q_post <- Q_spaz + Diagonal(N, 1 / sig2_err)

  chol_Q_post <- chol(Q_post)

  # Log-determinanti istantanei senza fattorizzazioni aggiuntive
  ldet_Q_spaz <- attr(Q_spaz, "log_det")
  ldet_Q_post <- 2 * sum(log(Matrix::diag(chol_Q_post)))

  # --- RISOLUZIONE VETTORIZZATA: ZERO CICLI FOR ---
  Z_mat <- Y_mat - curr_mu
  B_mat <- Z_mat / sig2_err

  Mu_post_mat <- as.matrix(solve(chol_Q_post, solve(t(chol_Q_post), B_mat)))

  err_sq_vec <- colSums((Z_mat - Mu_post_mat)^2) / sig2_err
  spat_energy_vec <- colSums(Mu_post_mat * (Q_spaz %*% Mu_post_mat))

  quad_sum <- sum(err_sq_vec + spat_energy_vec)

  log_lik <- -0.5 * N * M * log(sig2_err) + 0.5 * M * ldet_Q_spaz - 0.5 * M * ldet_Q_post - 0.5 * quad_sum
  return(log_lik)
}
