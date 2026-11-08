# ==============================================================================
# File: R/nngp_marginal.R
# Calcolo Stabile della Verosimiglianza Marginale (Identità di Rue & Held)
# ==============================================================================

log_marg_likelihood_robust <- function(sig2_spaz, phi, sig2_err, Y_mat, curr_mu, nngp_struct) {
  N <- nrow(Y_mat)
  M <- ncol(Y_mat)

  # 1. Matrici di Precisione
  Q_spaz <- build_nngp_precision(nngp_struct, sig2_spaz, phi)
  Q_post <- Q_spaz + Diagonal(N, 1 / sig2_err)

  # Cholesky della posteriore sparsa
  chol_Q_post <- chol(Q_post)

  # 2. Determinanti (log)
  ldet_Q_spaz <- as.numeric(determinant(Q_spaz, logarithm = TRUE)$modulus)
  ldet_Q_post <- as.numeric(determinant(Q_post, logarithm = TRUE)$modulus)

  # 3. Termini quadratici valutati alla media a posteriori (SENZA SOTTRAZIONI!)
  quad_sum <- 0
  for(i in 1:M) {
    z_i <- Y_mat[, i] - curr_mu
    b_i <- z_i / sig2_err

    # Calcoliamo la media a posteriori esatta per l'anno i
    mu_post_i <- as.vector(solve(chol_Q_post, solve(t(chol_Q_post), b_i)))

    # Termine A: Errore quadratico puro (Streptamente Positivo)
    err_sq <- sum((z_i - mu_post_i)^2) / sig2_err

    # Termine B: Energia spaziale (Streptamente Positivo)
    spat_energy <- as.numeric(t(mu_post_i) %*% Q_spaz %*% mu_post_i)

    # Li SOMMIAMO! Niente cancellazione catastrofica
    quad_sum <- quad_sum + (err_sq + spat_energy)
  }

  # 4. Assembliamo la Log-Verosimiglianza di Rue & Held
  log_lik <- -0.5 * N * M * log(sig2_err) + 0.5 * M * ldet_Q_spaz - 0.5 * M * ldet_Q_post - 0.5 * quad_sum
  return(log_lik)
}
