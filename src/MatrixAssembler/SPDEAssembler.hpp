#pragma once

#include <RcppEigen.h>


enum class SPDEType{
    ISOTROPIC = 0,
    ANISOTROPIC = 1
};

/**
 * @class SPDEAssembler
 * @brief Abstract base class for assembling Finite Element discretization matrices
 *        used in Stochastic Partial Differential Equation (SPDE) models.
 *
 * Stores mesh geometry, connectivity, and differential operator matrices.
 * @param V Dense vertex matrix of the mesh (N_v X d), where each row contains the spatial coordinates of a vertex.
 * @param F Dense face/element matrix of the mesh (N_f \times 3), where each row contains the indices of vertices in V that form a triangle. 
 * @property Kx Sparse stiffness matrix along the x-direction.
 * @property Ky Sparse stiffness matrix along the y-direction.
 * @property C  Sparse mass/capacity matrix.
 */
class SPDEAssembler {    
protected:
    Eigen::MatrixXd V;  ///< Vertex matrix of the mesh.
    Eigen::MatrixXi F;  ///< Faces matrix of the mesh. 

    Eigen::SparseMatrix<double> Kx; ///< Stiffness matrix (x-direction).
    Eigen::SparseMatrix<double> Ky; ///< Stiffness matrix (y-direction).
    Eigen::SparseMatrix<double> K; ///< Global Stifness Matrix
    Eigen::SparseMatrix<double> C;  ///< Mass / capacity matrix.


    SPDEType type;

public:
    /**
     * @brief Constructs an SPDEAssembler instance.
     * 
     * @param V_ Dense matrix containing vertex coordinates.
     * @param F_ Dense matrix containing edge connectivity indices.
     */
    SPDEAssembler(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_, SPDEType type_)
        : V(V_), F(F_), type(type_) {}

    virtual ~SPDEAssembler() = default;

    /**
     * @brief Assembles the differential and mass sparse matrices (Kx, Ky, C).
     */
    virtual void assemble_matrices() = 0;

    /// @brief Get the stiffness matrix along the x-direction. Empty if isotropic
    const Eigen::SparseMatrix<double>& get_Kx() const { return Kx; }
    /// @brief Get the stiffness matrix along the y-direction. Empty if isotropic
    const Eigen::SparseMatrix<double>& get_Ky() const { return Ky; }
    
    const Eigen::SparseMatrix<double>& get_K() const { return K; }
    const Eigen::SparseMatrix<double>& get_C()  const { return C; }
    SPDEType get_model_type() const { return type; }
};