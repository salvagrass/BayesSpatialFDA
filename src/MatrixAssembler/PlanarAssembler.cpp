#include "PlanarAssembler.hpp"
#include <igl/massmatrix.h>
#include<igl/cotmatrix.h>
#include<igl/grad.h>

void PlanarAssembler::assemble_matrices(){
    
    // Constructing the diagonal lumped mass matrix via libigl
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,C);
        

    // Constructing the stifness matrices
    
    // Computing the areas of the triangles
    Eigen::VectorXd double_areas;
    igl::doublearea(V,F,double_areas);

    Eigen::VectorXd areas = double_areas/2.0;

    // Computing the grad matrix, this matrix has 3F rows, three rows for each triangle and V columns, one for each vertex
    // For the first F rows G_ij holds the constant value of the x-direction of the gradient of the j-th basis function on i-th triangle
    // For the middle F rows G_ij holds the constant value of the y-direction of the gradient of the j-th basis function on i-th triangle
    Eigen::SparseMatrix<double> Gx,Gy,G;
    igl::grad(V,F,G);

    Gx = G.topRows(F.rows());
    Gy = G.middleRows(F.rows(),F.rows());

    // Our formula is K_{z,ij} = \sum_f^F \int_f grad(z,varphi_i) \cdot grad(z,varphi_j), since out gradients are constants along a triangle we can take them out
    // we have \sum_f^F grad(z,varphi_i) \grad(z,varphi_j) * Area(F) = \sum_f^F G_z(f,i)areas(f)G_z(f,j) = \sum_f^F G_z(f,i)diag(areas)_{f,f}G_z(f,j)
    // since the middle is a diagonal matrix we can write the sum as: \sum_w^F G_z(w,i)diag(areas)_(w,f) = (G_z^T*diag(areas))_(i,f) multiplying this
    // and summing on f gives us K_z = G_z^T*diag(areas)*G_z
    K11 = Gx.transpose()*areas.asDiagonal()*Gx;
    K22 = Gy.transpose()*areas.asDiagonal()*Gy;

    // As for the off diagonal stifness matrix the formula is K_{i,j} = \sum_f^F \int_f grad(varphi_i)_x grad(varphi_j)_y + grad(varphi_i)_y grad(varphi_j)_x
    // Using the same trick as before on a face of the triangle this are constants, so: 
    // K_{i,j} = \sum_f^F (grad(varphi_i)_x grad(varphi_j)_y + grad(varphi_i)_y*grad(varphi_j)_x)*Area(f)
    Eigen::SparseMatrix<double> K12_half = Gx.transpose()*areas.asDiagonal()*Gy;

    Eigen::SparseMatrix<double> K12_half_T = K12_half.transpose();

    K12 = K12_half + K12_half_T;
}


