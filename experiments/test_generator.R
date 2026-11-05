# File: experiments/test_generate_toy_data.R
devtools::load_all()

# 1. Mesh veloce giusto per avere i nodi
set.seed(42)
locs <- matrix(runif(200), ncol = 2)
mesh <- fm_mesh_2d(loc = locs, max.edge = c(0.1, 0.5))
# plot(mesh)
N <- nrow(mesh$loc)
M <- 10 # Numero di anni da simulare

# 2. Chiamata alla funzione con parametri noti
sim <- generate_toy_data(
  mesh = mesh,
  num_year = M,
  a_prior_mean = 5.0, # Ci aspettiamo media ~ 5.0 su M anni
  a_prior_sd = 1.0,
  mu_prior_mean = 10.0, # Ci aspettiamo ~ 10.0 sui nodi
  Q_mu = 1, # Precisione dummy per mu
  sigma_e = 0.5 # Varianza attesa = 0.25
)

# 3. Check visivo veloce dei risultati
cat("=== SANITY CHECK GENERATORE MULTI-ANNO ===\n\n")

cat("0. Dimensioni delle matrici (Y_mat, M x N):\n")
cat(
  "   Attese:",
  N,
  "x",
  M,
  "| Ottenute:",
  nrow(sim$y_mat),
  "x",
  ncol(sim$y_mat),
  "\n\n"
)

cat("1. Effetto temporale medio (a_true su", M, "anni):\n")
cat("   Atteso: ~ 5.0 | Ottenuto:", mean(sim$a_true), "\n\n")

cat("2. Media del campo spaziale globale (mu_true):\n")
cat("   Attesa: ~ 10.0 | Ottenuta:", mean(sim$mu_true), "\n\n")

cat("3. Media dell'anomalia globale (beta_tilde_mat):\n")
cat("   Attesa: ~ 0.0 | Ottenuta:", mean(sim$beta_tilde_mat), "\n\n")

cat("4. Varianza del rumore di misurazione (epsilon):\n")
# Vettorizziamo le matrici per calcolare la varianza globale dell'errore
var_epsilon_empirica <- var(as.vector(sim$y_mat - sim$latent_field_mat))
cat("   Attesa: ~ 0.25 | Ottenuta:", var_epsilon_empirica, "\n\n")
