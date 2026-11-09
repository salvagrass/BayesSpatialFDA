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

  ldet_Q_spaz <- as.numeric(determinant(Q_spaz, logarithm = TRUE)$modulus)
  ldet_Q_post <- as.numeric(determinant(Q_post, logarithm = TRUE)$modulus)

  # --- RISOLUZIONE VETTORIZZATA: ZERO CICLI FOR ---
  # Z_mat e B_mat coprono tutto il dataset
  Z_mat <- Y_mat - curr_mu
  B_mat <- Z_mat / sig2_err

  # solve() su un oggetto dCHMsimpl di Matrix accetta RHS multipli (matrici).
  # Questo esegue il C-level solve in blocco su tutto il campione simultaneamente!
  Mu_post_mat <- as.matrix(solve(chol_Q_post, solve(t(chol_Q_post), B_mat)))

  # Termine A: colSums di matrici evita i loop iterativi sugli M periodi
  err_sq_vec <- colSums((Z_mat - Mu_post_mat)^2) / sig2_err

  # Termine B: calcolo dell'energia spaziale vettorizzata
  # Invece di fare M prodotti t(x) %*% Q %*% x, usiamo il prodotto Hadamard e colSums
  spat_energy_vec <- colSums(Mu_post_mat * (Q_spaz %*% Mu_post_mat))

  # Sommiamo tutti i termini M generati (nessuna cancellazione numerica)
  quad_sum <- sum(err_sq_vec + spat_energy_vec)

  log_lik <- -0.5 * N * M * log(sig2_err) + 0.5 * M * ldet_Q_spaz - 0.5 * M * ldet_Q_post - 0.5 * quad_sum
  return(log_lik)
}
