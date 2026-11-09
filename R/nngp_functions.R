# ==============================================================================
# MODULO 1: COSTRUZIONE DEL GRAFO NNGP (Nearest Neighbors) OTTIMIZZATO
# ==============================================================================

#' Costruisce la struttura dei vicini e PRE-CALCOLA tutte le distanze interne
build_nngp_neighbors <- function(coords, m = 15) {
  N <- nrow(coords)

  # 1. Ordinamento spaziale
  ord <- order(coords[, 1], coords[, 2])
  coords_ord <- coords[ord, ]

  nn_idx <- matrix(NA, nrow = N, ncol = m)
  nn_dist <- matrix(NA, nrow = N, ncol = m)

  # NOVITÀ: Lista per conservare in cache le matrici di distanza INTERNE ai vicini
  # Questo evita di chiamare dist() migliaia di volte durante la MCMC!
  nn_dist_matrix_list <- vector("list", N)

  cat("Costruzione struttura NNGP e Caching per", N, "vertici (m =", m, ")...\n")

  # 2. Ricerca e Caching
  for (i in 2:N) {
    num_neighbors <- min(m, i - 1)
    prev_coords <- coords_ord[1:(i - 1), , drop = FALSE]

    # Distanza dal punto i ai precedenti
    dist_to_prev <- sqrt((prev_coords[, 1] - coords_ord[i, 1])^2 +
                           (prev_coords[, 2] - coords_ord[i, 2])^2)

    nearest <- order(dist_to_prev)[1:num_neighbors]

    # Salviamo indici e distanze (dal punto i ai vicini)
    nn_idx[i, 1:num_neighbors] <- nearest
    nn_dist[i, 1:num_neighbors] <- dist_to_prev[nearest]

    # CACHING: Calcoliamo QUI (una volta sola per tutta la vita del modello)
    # le distanze reciproche TRA i vicini.
    if(num_neighbors > 1) {
      coords_neigh <- coords_ord[nearest, , drop = FALSE]
      nn_dist_matrix_list[[i]] <- as.matrix(dist(coords_neigh))
    } else {
      nn_dist_matrix_list[[i]] <- matrix(0, 1, 1) # Nodo 2 ha 1 solo vicino, distanza da se stesso = 0
    }
  }

  return(list(
    coords_ord = coords_ord,
    ord_indices = ord,
    nn_idx = nn_idx,
    nn_dist = nn_dist,
    nn_dist_matrix = nn_dist_matrix_list, # <--- La chiave della nuova velocità
    m = m
  ))
}

# ==============================================================================
# MODULO 2: COSTRUZIONE DELLA MATRICE DI PRECISIONE NNGP (ULTRA VELOCE)
# ==============================================================================

#' Costruisce la matrice di precisione sfruttando il caching delle distanze
build_nngp_precision <- function(nngp_struct, sigma2, phi) {

  N <- nrow(nngp_struct$coords_ord)
  m <- nngp_struct$m
  nn_idx <- nngp_struct$nn_idx

  max_elements <- N + N * m
  i_vec <- numeric(max_elements)
  j_vec <- numeric(max_elements)
  x_vec <- numeric(max_elements)

  F_vec <- numeric(N)

  # Nodo 1
  F_vec[1] <- sigma2
  i_vec[1] <- 1
  j_vec[1] <- 1
  x_vec[1] <- 1
  idx_count <- 2 # Inizializziamo a 2 per il prossimo inserimento

  # Ciclo sui restanti nodi (Senza chiamate a dist() e senza cicli for annidati)
  for (i in 2:N) {
    num_neigh <- sum(!is.na(nn_idx[i, ]))
    neigh <- nn_idx[i, 1:num_neigh]

    # --- Calcolo Covarianze ULTRA-VELOCE tramite Cache ---
    # Non usiamo più dist(coords_neigh). Leggiamo dalla lista!
    dist_neigh <- nngp_struct$nn_dist_matrix[[i]]
    C_neigh <- sigma2 * exp(-phi * dist_neigh)

    dist_i_neigh <- nngp_struct$nn_dist[i, 1:num_neigh]
    c_cross <- sigma2 * exp(-phi * dist_i_neigh)

    # --- Risoluzione micro-sistema lineare ---
    B_i <- as.vector(solve(C_neigh, c_cross))
    F_vec[i] <- sigma2 - sum(c_cross * B_i)

    # --- Assegnazione Vettorizzata nella matrice sparsa (NESSUN CICLO FOR) ---
    end_idx <- idx_count + num_neigh

    # 1. Mettiamo l'1 sulla diagonale
    i_vec[idx_count] <- i
    j_vec[idx_count] <- i
    x_vec[idx_count] <- 1

    # 2. Mettiamo i coefficienti -B_i nei posti giusti IN BLOCCO
    i_vec[(idx_count + 1):end_idx] <- i
    j_vec[(idx_count + 1):end_idx] <- neigh
    x_vec[(idx_count + 1):end_idx] <- -B_i

    # Aggiorniamo il contatore globale
    idx_count <- end_idx + 1
  }

  valid_len <- idx_count - 1
  i_vec <- i_vec[1:valid_len]
  j_vec <- j_vec[1:valid_len]
  x_vec <- x_vec[1:valid_len]

  B_nngp <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(N, N))
  F_inv <- Matrix::Diagonal(x = 1 / F_vec)

  Q_nngp <- Matrix::t(B_nngp) %*% F_inv %*% B_nngp

  # Ordinamento finale
  ord <- nngp_struct$ord_indices
  inv_ord <- order(ord)
  return(Q_nngp[inv_ord, inv_ord])
}
