#pragma once

#include <RcppEigen.h>


/**
 * @class SPDEAssembler
 * @brief Abstract base class for assembling Finite Element discretization matrices
 *        used in Stochastic Partial Differential Equation (SPDE) models.
 *
 * Stores mesh geometry, connectivity, and differential operator matrices.
 * @param V Dense vertex matrix of the mesh (N_v X d), where each row contains the spatial coordinates of a vertex.
 * @param F Dense face/element matrix of the mesh (N_f \times 3), where each row contains the indices of vertices in V that form a triangle. 
 * @property K11 Sparse stiffness matrix along the x-direction.
 * @property K22 Sparse stiffness matrix along the y-direction.
 * @property K12 Cross diffusion matrix
 * @property C  Sparse mass/capacity matrix.
 */
class SPDEAssembler {    
protected:
    Eigen::MatrixXd V;  ///< Vertex matrix of the mesh.
    Eigen::MatrixXi F;  ///< Faces matrix of the mesh. 

    Eigen::SparseMatrix<double> K11; ///< Stiffness matrix (x-direction).
    Eigen::SparseMatrix<double> K22; ///< Stiffness matrix (y-direction).
    Eigen::SparseMatrix<double> K12; ///< Cross diffusion matrix / Cross Coupling stifness matrix

    Eigen::SparseMatrix<double> C;  ///< Mass / capacity matrix

public:
    /**
     * @brief Constructs an SPDEAssembler instance.
     * 
     * @param V_ Dense matrix containing vertex coordinates.
     * @param F_ Dense matrix containing edge connectivity indices.
     */
    SPDEAssembler(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_)
        : V(V_), F(F_) {}

    virtual ~SPDEAssembler() = default;

    /**
     * @brief Assembles the differential and mass sparse matrices (Kx, Ky, C).
     */
    virtual void assemble_matrices() = 0;

    /// @brief Get the stiffness matrix along the x-direction.
    const Eigen::SparseMatrix<double>& get_K11() const { return K11; }
    /// @brief Get the stiffness matrix along the y-direction. 
    const Eigen::SparseMatrix<double>& get_K22() const { return K22; }
    /// @brief Get the cross product stifness matrix, associated with the off diagonal term
    const Eigen::SparseMatrix<double>& get_K12() const {return K12; }
    
    const Eigen::SparseMatrix<double>& get_C()  const { return C; }
};