# ==============================================================================
# MODULO 1: COSTRUZIONE DEL GRAFO NNGP (Nearest Neighbors)
# ==============================================================================

#' Costruisce la struttura dei vicini per il modello NNGP
#'
#' @param coords Matrice Nx2 delle coordinate spaziali
#' @param m Numero massimo di vicini da considerare (solitamente 10-15)
#' @return Una lista contenente le coordinate ordinate, gli indici di ordinamento,
#'         la matrice degli indici dei vicini e la matrice delle distanze.
build_nngp_neighbors <- function(coords, m = 15) {
  N <- nrow(coords)

  # 1. Ordinamento spaziale (Sweep lungo l'asse X, poi Y in caso di parità)
  # Fondamentale per garantire che la matrice di precisione Q sia triangolare
  ord <- order(coords[, 1], coords[, 2])
  coords_ord <- coords[ord, ]

  # Strutture dati: matrici N x m.
  # La riga i conterrà gli indici (e le distanze) dei vicini del punto i.
  nn_idx <- matrix(NA, nrow = N, ncol = m)
  nn_dist <- matrix(NA, nrow = N, ncol = m)

  cat("Costruzione struttura NNGP per", N, "vertici con m =", m, "vicini...\n")

  # 2. Ricerca dei vicini (Condizionamento all'indietro)
  # Il punto 1 non ha vicini precedenti. Iniziamo dal punto 2.
  for (i in 2:N) {
    # Quanti vicini precedenti esistono al massimo per questo punto?
    num_neighbors <- min(m, i - 1)

    # Estraiamo le coordinate di tutti i punti precedenti
    prev_coords <- coords_ord[1:(i - 1), , drop = FALSE]

    # Calcoliamo la distanza euclidea tra il punto i e tutti i punti precedenti
    # (Vettorizzato per velocità in base R)
    dist_to_prev <- sqrt((prev_coords[, 1] - coords_ord[i, 1])^2 +
                           (prev_coords[, 2] - coords_ord[i, 2])^2)

    # Ordiniamo le distanze per trovare i più vicini
    nearest <- order(dist_to_prev)[1:num_neighbors]

    # Salviamo i risultati (gli indici si riferiscono all'array ordinato)
    nn_idx[i, 1:num_neighbors] <- nearest
    nn_dist[i, 1:num_neighbors] <- dist_to_prev[nearest]
  }

  return(list(
    coords_ord = coords_ord,
    ord_indices = ord,
    nn_idx = nn_idx,     # Chi sono i miei vicini?
    nn_dist = nn_dist,   # Quanto distano da me?
    m = m
  ))
}

# ==============================================================================
# MODULO 2: COSTRUZIONE DELLA MATRICE DI PRECISIONE NNGP (SPARSA)
# ==============================================================================

#' Costruisce la matrice di precisione sparsa Q usando il grafo NNGP
#'
#' @param nngp_struct Lista restituita da build_nngp_neighbors()
#' @param sigma2 Varianza spaziale (marginal variance)
#' @param phi Parametro di decadimento spaziale (spatial decay, range = 3/phi)
#' @return Matrice di precisione sparsa Q_final nelle coordinate originali
build_nngp_precision <- function(nngp_struct, sigma2, phi) {

  N <- nrow(nngp_struct$coords_ord)
  m <- nngp_struct$m
  nn_idx <- nngp_struct$nn_idx
  coords <- nngp_struct$coords_ord

  # Strutture per creare la matrice sparsa in formato "triplet" (i, j, valore)
  # B_nngp ha 1 sulla diagonale e fino a m valori fuori diagonale per riga
  max_elements <- N + N * m
  i_vec <- numeric(max_elements)
  j_vec <- numeric(max_elements)
  x_vec <- numeric(max_elements)

  F_vec <- numeric(N)
  idx_count <- 1

  # Nodo 1: non ha vicini precedenti, la sua varianza condizionata è quella marginale
  F_vec[1] <- sigma2
  i_vec[idx_count] <- 1
  j_vec[idx_count] <- 1
  x_vec[idx_count] <- 1
  idx_count <- idx_count + 1

  # Ciclo sui restanti nodi
  for (i in 2:N) {
    # Identifichiamo i vicini effettivi del nodo (ignorando gli NA iniziali)
    num_neigh <- sum(!is.na(nn_idx[i, ]))
    neigh <- nn_idx[i, 1:num_neigh]

    # Matrice di covarianza tra i vicini (C_neigh)
    # È qui che risparmiamo: questa matrice è al massimo m x m (es. 15x15)!
    coords_neigh <- coords[neigh, , drop = FALSE]
    dist_neigh <- as.matrix(dist(coords_neigh))
    C_neigh <- sigma2 * exp(-phi * dist_neigh)

    # Covarianza incrociata tra il nodo i e i suoi vicini (c_cross)
    # nngp_struct$nn_dist contiene già le distanze calcolate nello step precedente
    dist_i_neigh <- nngp_struct$nn_dist[i, 1:num_neigh]
    c_cross <- sigma2 * exp(-phi * dist_i_neigh)

    # Risolviamo il micro-sistema lineare C_neigh * B_i = c_cross
    # solve() su una 15x15 è istantaneo
    B_i <- as.vector(solve(C_neigh, c_cross))

    # Calcoliamo la varianza residua (condizionata)
    F_vec[i] <- sigma2 - sum(c_cross * B_i)

    # --- Popoliamo la matrice sparsa B_nngp ---
    # Diagonale = 1
    i_vec[idx_count] <- i
    j_vec[idx_count] <- i
    x_vec[idx_count] <- 1
    idx_count <- idx_count + 1

    # Fuori diagonale (i pesi dei vicini col segno meno)
    for(k in 1:num_neigh) {
      i_vec[idx_count] <- i
      j_vec[idx_count] <- neigh[k]
      x_vec[idx_count] <- -B_i[k]
      idx_count <- idx_count + 1
    }
  }

  # Tronchiamo i vettori alla loro lunghezza effettiva
  valid_len <- idx_count - 1
  i_vec <- i_vec[1:valid_len]
  j_vec <- j_vec[1:valid_len]
  x_vec <- x_vec[1:valid_len]

  # Creiamo la matrice sparsa B usando il pacchetto Matrix
  B_nngp <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(N, N))

  # Creiamo la matrice diagonale delle precisioni condizionate
  F_inv <- Matrix::Diagonal(x = 1 / F_vec)

  # Ricaviamo la precisione globale (ordinata)
  Q_nngp <- Matrix::t(B_nngp) %*% F_inv %*% B_nngp

  # ULTIMO STEP: Riportiamo la matrice Q all'ordinamento originale dei dati
  ord <- nngp_struct$ord_indices
  inv_ord <- order(ord)
  Q_final <- Q_nngp[inv_ord, inv_ord]

  return(Q_final)
}
