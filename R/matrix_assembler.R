#' Build SPDE Matrices (Isotropic)
#'
#' This function takes a spatial mesh and constructs the finite element 
#' matrices (Mass matrix C and Stiffness matrix K) required for SPDE modeling.
#'
#' @param mesh An `fm_mesh_2d` object created via `fmesher`.
#' @param manifold Character string. Either `"planar"` or `"sphere"`.
#'
#' @return A list containing two sparse matrices: `C` (Mass) and `K` (Stiffness).
#' @export
#'
#' @examples
#' # mesh <- fmesher::fm_rcdt_2d_inla(...)
#' # matrices <- build_spde_matrices(mesh, manifold = "planar")
build_spde_matrices <- function(mesh, manifold = c("planar", "sphere")) {
  
  manifold <- match.arg(manifold)
  
  if (!inherits(mesh, "fm_mesh_2d")) {
    stop("Input 'mesh' must be an 'fm_mesh_2d' object from the fmesher package.")
  }
  
  # 2. Mappatura Enum per il C++
  manifold_int <- if (manifold == "planar") 0L else 1L
  
  # 3. Estrazione dati puliti dalla mesh complessa
  V <- mesh$loc          # Matrice delle coordinate (Nx2 o Nx3)
  F_r <- mesh$graph$tv   # Matrice dei triangoli (1-based index)
  
  # 4. IL PASSAGGIO CRITICO: R conta da 1, C++ conta da 0
  F_cpp <- F_r - 1L
  
  # 5. Chiamata al motore C++ invisibile
  # Usiamo .Call o la funzione mappata da Rcpp
  matrices <- .build_spde_matrices_cpp(V, F_cpp)
  
  # 6. Ritorno pulito
  return(matrices)
}