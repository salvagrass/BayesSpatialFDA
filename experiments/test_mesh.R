devtools::load_all()
library(fmesher)

# Creating a lattice
lattice <- fmesher::fm_lattice_2d(x = c(0, 1), y = c(0, 1))

# Mesh Generation
mesh_crossed <- fmesher::fm_rcdt_2d_inla(lattice = lattice, extend = FALSE)

plot(mesh_crossed, main = "Cross Structured Mesh")

# Testing C++ function
matrici_test <- build_spde_matrices(mesh_crossed, manifold = "planar")
print(matrici_test$K)
print(matrici_test$C)

# Computing same matrices with fmeshere
fem_fmesher <- fmesher::fm_fem(mesh_crossed)


# Computing the difference between the two matrices
max(abs(matrici_test$C - fem_fmesher$c0))
max(abs(matrici_test$K - fem_fmesher$g1))
