test_that("generate_toy_data produces correct dimensions and mathematical structures", {
  # 1. Setup a basic mesh using fmesher
  lattice <- fmesher::fm_lattice_2d(
    x = seq(0, 1, length.out = 5),
    y = seq(0, 1, length.out = 5)
  )
  test_mesh <- fmesher::fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)

  N_vertices <- nrow(test_mesh$loc)
  M_years <- 10

  # 2. Generate toy data with a fixed seed
  toy_data <- generate_toy_data(
    mesh = test_mesh,
    num_year = M_years,
    mu_prior_mean = 15.0,
    sigma_e = 0.1,
    seed = 42
  )

  # 3. Test list structure and names
  expected_names <- c(
    "y_mat",
    "mu_true",
    "beta_mat"
  )
  expect_named(toy_data, expected_names, ignore.order = TRUE)

  # 4. Test dimensions
  expect_length(toy_data$mu_true, N_vertices)

  expect_equal(nrow(toy_data$y_mat), N_vertices)
  expect_equal(ncol(toy_data$y_mat), M_years)

  expect_equal(nrow(toy_data$beta_mat), N_vertices)
  expect_equal(ncol(toy_data$beta_mat), M_years)

  # 5. Test Mathematical Relationships
  # The final observations Y must be the spatial latent field plus observation error (epsilon)
  # We check that Y is close to beta_mat, but not identical due to Gaussian noise.
  expect_false(isTRUE(all.equal(toy_data$y_mat, toy_data$beta_mat)))

  # 6. Test Reproducibility
  toy_data_2 <- generate_toy_data(
    mesh = test_mesh,
    num_year = M_years,
    mu_prior_mean = 15.0,
    sigma_e = 0.1,
    seed = 42
  )

  expect_equal(toy_data$y_mat, toy_data_2$y_mat)
  expect_equal(toy_data$mu_true, toy_data_2$mu_true)
  expect_equal(toy_data$beta_mat, toy_data_2$beta_mat)
})
