library(fmesher)
library(Matrix)
library(ggplot2)
library(dplyr)

# 1. Setting up the mesh
lattice <- fm_lattice_2d(
  x = seq(0, 10, length.out = 30),
  y = seq(0, 10, length.out = 30)
)
mesh <- fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)
N <- nrow(mesh$loc)

# 2. Finding the central node of the mesh
center_x <- 5
center_y <- 5
distances <- sqrt((mesh$loc[, 1] - center_x)^2 + (mesh$loc[, 2] - center_y)^2)
central_node <- which.min(distances)
central_df <- data.frame(
  x = mesh$loc[central_node, 1],
  y = mesh$loc[central_node, 2]
)

# 3. Defining scenarios for anisotropy
scenarios <- list(
  list(gamma = 1.0, theta = 0.0, title = "1. Isotropic (gamma = 1.0)"),
  list(
    gamma = 4.0,
    theta = 0.0,
    title = "2. Anisotropic: Stretched X (gamma=4)"
  ),
  list(
    gamma = 4.0,
    theta = pi / 4,
    title = "3. Anisotropic: Rotated 45 degrees"
  ),
  list(
    gamma = 0.25,
    theta = 0.0,
    title = "4. Anisotropic: Stretched Y (gamma=0.25)"
  )
)

triangles <- mesh$graph$tv

# 4. Computing emperical covariance in the field
plot_data_list <- list()

for (scen in scenarios) {
  data <- generate_toy_data(
    mesh,
    num_year = 300,
    gamma = scen$gamma,
    theta = scen$theta,
    seed = 42
  )

  central_series <- data$beta_mat[central_node, ]
  emp_cov <- apply(data$beta_mat, 1, function(node_series) {
    cor(node_series, central_series)
  })

  emp_cov[is.na(emp_cov)] <- 0
  tri_values <- rowMeans(matrix(emp_cov[triangles], ncol = 3))

  for (j in seq_len(nrow(triangles))) {
    node_ids <- triangles[j, ]
    tri_df <- data.frame(
      x = mesh$loc[node_ids, 1],
      y = mesh$loc[node_ids, 2],
      triangle_id = j,
      value = tri_values[j],
      scenario = scen$title
    )
    plot_data_list[[length(plot_data_list) + 1]] <- tri_df
  }
}

plot_data <- bind_rows(plot_data_list)
plot_data$scenario <- factor(
  plot_data$scenario,
  levels = sapply(scenarios, function(s) s$title)
)

# 5. Plotting Empirical Covariance
p <- ggplot(plot_data) +
  geom_polygon(
    aes(x = x, y = y, group = interaction(scenario, triangle_id), fill = value),
    color = NA
  ) +
  geom_point(data = central_df, aes(x = x, y = y), color = "black", size = 2) +
  scale_fill_viridis_c(option = "plasma") + # Scala automatica basata sui dati reali
  coord_fixed(ratio = 1) +
  facet_wrap(~scenario, nrow = 2, labeller = label_wrap_gen(width = 30)) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Empirical Spatial Correlation from Central Reference Node",
    subtitle = "Validating SPDE Anisotropy Formulation (Gamma & Theta)",
    fill = "Correlation",
    x = "X Coordinate",
    y = "Y Coordinate"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

print(p)
