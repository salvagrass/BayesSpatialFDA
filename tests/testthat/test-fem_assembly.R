test_that("Testing assembly of the mesh using libigl, comparing with fmesher", {
  # 1. Creating a complex lattice (10x10 grid = 100 vertices, 162 triangles)
  lattice <- fmesher::fm_lattice_2d(
    x = seq(0, 1, length.out = 10),
    y = seq(0, 1, length.out = 10)
  )

  # 2. Mesh generation
  mesh_crossed <- fmesher::fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)

  # 3. Testing C++ function (returns C, K11, K22, and K12)
  test_matrixes <- build_spde_matrices(mesh_crossed, manifold = "planar")

  # 4. Computing reference matrices with fmesher
  fem_fmesher <- fmesher::fm_fem(mesh_crossed)

  # 5. Cleaning f_mesher matrices and converting to dgCMatrix
  C_fmesher <- as(fem_fmesher$c0, "CsparseMatrix")
  K_fmesher <- as(fem_fmesher$g1, "CsparseMatrix")

  # --- TEST 1: STANDARD MATRICES MATCH FMESHER ---

  # Check if Lumped Mass Matrix (C) matches
  expect_equal(
    test_matrixes$C,
    C_fmesher,
    tolerance = 1e-10,
    info = "Lumped mass matrix (C) should match fmesher exactly."
  )

  # Check if Isotropic Stiffness Matrix (K11 + K22) matches fmesher's g1
  K_iso_test <- test_matrixes$K11 + test_matrixes$K22
  expect_equal(
    K_iso_test,
    K_fmesher,
    tolerance = 1e-10,
    info = "The sum K11 + K22 must exactly equal the isotropic stiffness matrix (g1) from fmesher."
  )

  # --- TEST 2: ALGEBRAIC INVARIANTS FOR K12 (CROSS-DERIVATIVE MATRIX) ---

  # A. Symmetry Check
  # K12 = Gx^T * A * Gy + Gy^T * A * Gx must mathematically be symmetric.
  expect_true(
    isSymmetric(test_matrixes$K12),
    info = "K12 (cross-derivative stiffness matrix) must be symmetric."
  )

  # B. Constant Field Derivative Check (Zero Row/Col Sums)
  # The derivative of a constant field is zero.
  # Therefore, the sum of each row and each column must equal zero (within floating-point precision).
  row_sums_K12 <- rowSums(test_matrixes$K12)
  expect_true(
    all(abs(row_sums_K12) < 1e-12),
    info = "Row sums of K12 must be strictly zero."
  )
})
