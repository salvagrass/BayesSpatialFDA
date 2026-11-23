# Se stai sviluppando il pacchetto con devtools:
 devtools::load_all()
# Altrimenti:
# source("R/data_prep.R")

# 1) Caricamento e pulizia globale (controlla che esca esattamente N = 64442 x 25)
clean <- load_clean_climate_data("data/climate_data.RData")
dim(clean$y_mat)
# Atteso: [1] 64442    25

# 2) Test A: Quadrato 2000 x 2000 km centrato su Roma (con riscalamento unitario)
roma_sq <- extract_planar_region(
  clean,
  center = c(41.9028, 12.4964),
  shape = "square",
  side_km = 2000,
  scale_to_unit = TRUE
)
mesh_sq <- build_planar_mesh(roma_sq, max_edge = c(0.08, 0.25), offset = 0.2)

cat("Quadrato Roma -> Osservazioni N:", nrow(roma_sq$coords),
    "| Nodi mesh V:", mesh_sq$mesh$n,
    "| Dim Phi:", paste(dim(mesh_sq$Phi), collapse = "x"), "\n")

# Verifica che ogni osservazione cada dentro la mesh (somma righe di Phi == 1)
stopifnot(max(abs(Matrix::rowSums(mesh_sq$Phi) - 1)) < 1e-10)

# Verifica integrazione diretta con il vostro motore C++:
spde_mats <- build_spde_matrices(mesh_sq$mesh, manifold = "planar")
str(spde_mats)

# 3) Test B: KNN con K = 400 punti più vicini a Roma
roma_knn <- prepare_real_planar_data(
  path = "data/climate_data.RData",
  center = c(41.9028, 12.4964),
  shape = "knn",
  k = 400L,
  scale_to_unit = FALSE,
  offset = 300
)

# 4) Plot diagnostico affiancato (Quadrato vs KNN)
par(mfrow = c(1, 2))

plot(mesh_sq$mesh, main = sprintf("Quadrato Roma (N=%d, V=%d)",
                                  nrow(roma_sq$coords), mesh_sq$mesh$n))
polygon(roma_sq$boundary, border = "red", lwd = 2)
points(roma_sq$coords, pch = 16, cex = 0.5, col = "steelblue")

plot(roma_knn$mesh, main = sprintf("KNN Roma (K=%d, V=%d)",
                                   nrow(roma_knn$coords), roma_knn$mesh$n))
polygon(roma_knn$boundary, border = "red", lwd = 2)
points(roma_knn$coords, pch = 16, cex = 0.5, col = "darkgreen")
