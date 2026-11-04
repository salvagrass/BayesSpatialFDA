#include "PlanarAssembler.hpp"
#include <igl/massmatrix.h>
#include<igl/cotmatrix.h>
#include <vector>

void PlanarAssembler::assemble_matrices(){
    
    // Constructing the diagonal lumped mass matrix via libigl
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,C);
        
    //  Constructing the isotropic stifness matrix and align with sign convention in FEM
    igl::cotmatrix(V,F,K);
    K.coeffs() *= -1.0;
}