// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <memory>

#include "MatrixAssembler/SPDEAssembler.hpp"
#include "MatrixAssembler/PlanarAssembler.hpp"
// #include "MatrixAssembler/SphereAssembler.hpp" // Uncomment when ready

// [[Rcpp::export(name = ".build_spde_matrices_cpp")]]
Rcpp::List build_spde_matrices_cpp(const Eigen::MatrixXd &V,
                                   const Eigen::MatrixXi &F)
{
    // Smart pointer for zero-leak memory management
    std::unique_ptr<SPDEAssembler> assembler;

    // Will have a factory pattern to instantiate the correct geometry
    assembler = std::make_unique<PlanarAssembler>(V, F);
    // Execute the FEM heavy lifting
    assembler->assemble_matrices();

    // Return the list of sparse matrices (C and K for isotropic) to R
    return Rcpp::List::create(
        Rcpp::Named("C") = assembler->get_C(),
        Rcpp::Named("K11") = assembler->get_K11(),
        Rcpp::Named("K12") = assembler->get_K12(),
        Rcpp::Named("K22") = assembler->get_K22());
}