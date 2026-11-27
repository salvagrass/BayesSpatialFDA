/**
 * @file sampler_types.hpp
 * @brief Core data structures, preallocated MCMC workspace, and SPDE precision matrix assembly.
 *
 * Defines the state, read-only sufficient statistics, thread-local solver memory, and output containers.
 */
#pragma once

#include <cmath>
<<<<<<< HEAD
#include <cstdint>
=======
>>>>>>> main
#include <map>
#include <memory>
#include <random>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>


/**
 * @enum ThetaIdx
 * @brief Positional indices for the SPDE hyperparameter vector \f$ \theta = (\kappa, \sigma_u, v_1, v_2) \f$.
 */
enum ThetaIdx : int {
    KAPPA   = 0,
    SIGMA_BETA = 1,
    V1      = 2,
    V2      = 3
};

/**
 * @struct ParallelConfig
 * @brief Threading configuration across MCMC chains and linear algebra routines.
 */
struct ParallelConfig{
    int num_chains = 1;
    int chol_threads = 1;
<<<<<<< HEAD
    int year_threads = 1;
=======
>>>>>>> main
};


/**
 * @struct ModelState
 * @brief Mutable state of a single Markov chain at the current Gibbs iteration.
 */
struct  ModelState{
    Eigen::MatrixXd spatial_coeffs;
    Eigen::VectorXd spde_hyperparameters;
    Eigen::VectorXd gp_mean;
    double error_variance;
};



/**
 * @struct SamplerContext
 * @brief Immutable, read-only context shared safely across all parallel MCMC chains.
 *
 * Stores precomputed sufficient statistics (eliminating \f$ \mathcal{O}(N) \f$ grid operations inside
 * the MCMC loop), finite element geometric matrices on \f$ \mathbb{S}^2 \f$, and problem dimensions.
 */
struct SamplerContext{
<<<<<<< HEAD
    const Eigen::SparseMatrix<double> PhiTPhi;
    const Eigen::MatrixXd PhiTY;
    const Eigen::VectorXd PhiTY_sum; //rowise sum of PhiTY

=======

    SamplerContext(
        Eigen::SparseMatrix<double> PhiTPhi_,
        Eigen::MatrixXd PhiTY_,
        Eigen::VectorXd C_,
        Eigen::SparseMatrix<double> Gx_,
        Eigen::SparseMatrix<double> Gy_,
        Eigen::SparseMatrix<double> Gxy_,
        const double YTY_sum_,
        const size_t num_iters_,
        const ParallelConfig config_,
        const long num_observations_
    ) : PhiTPhi(std::move(PhiTPhi_)), 
    PhiTY(std::move(PhiTY_)),
    PhiTY_sum(PhiTY.rowwise().sum()), 
    C(std::move(C_)),
    C_inv(C.cwiseInverse()),
    Gx(std::move(Gx_)), 
    Gy(std::move(Gy_)),
    Gxy(std::move(Gxy_)), 
    YTY_sum(YTY_sum_), 
    num_observations(num_observations_),
    num_vertices(PhiTPhi.rows()),
    num_years(PhiTY.cols()),
    Identity_V([this](){
            Eigen::SparseMatrix<double> id(num_vertices,num_vertices);
            id.setIdentity();
            return id;
        }()),
    num_iters(num_iters_),
    config(config_){}


    const Eigen::SparseMatrix<double> PhiTPhi;
    const Eigen::MatrixXd PhiTY;
    const Eigen::VectorXd PhiTY_sum; //rowise sum of PhiTY
    
>>>>>>> main

    const Eigen::VectorXd C; // Since C is diagonal we store it as a vector
    const Eigen::VectorXd C_inv; // We also store the inverse of C

    const Eigen::SparseMatrix<double> Gx;
    const Eigen::SparseMatrix<double> Gy;
    const Eigen::SparseMatrix<double> Gxy;

<<<<<<< HEAD
    const double YTY;
    
    const int num_observations = 0;
    const int num_vertices = 0;
    const int num_years = 0;
=======
    const double YTY_sum;
    
    const Eigen::Index num_observations;
    const Eigen::Index num_vertices;
    const Eigen::Index num_years;

    const Eigen::SparseMatrix<double> Identity_V; 

    const size_t num_iters;

    const ParallelConfig config;
>>>>>>> main
};

/**
 * @struct SamplerMemory
 * @brief Preallocated, thread-local workspace for a single MCMC chain.
 *
 * Eliminates dynamic heap allocations during Gibbs sweeps by holding reusable sparse
 * buffers, dense scratchpads, and symbolically analyzed sparse Cholesky solvers grouped
 * into 2-ring (SPDE / posterior latent field) and 1-ring (non-centered ASIS) sparsity patterns.
 */
struct SamplerMemory {
    using SpMat  = Eigen::SparseMatrix<double>;
    using Solver = Eigen::SimplicialLLT<SpMat>;

<<<<<<< HEAD
    SpMat Q_spde_curr, Q_spde_prop;
    SpMat buffer_1ring, buffer_2ring;

    std::unique_ptr<Solver> mu_solver   = std::make_unique<Solver>();
    std::unique_ptr<Solver> asis_solver = std::make_unique<Solver>();
    std::unique_ptr<Solver> mh_curr     = std::make_unique<Solver>();
    std::unique_ptr<Solver> mh_prop     = std::make_unique<Solver>();
    std::unique_ptr<Solver> beta_curr   = std::make_unique<Solver>();
    std::unique_ptr<Solver> beta_prop   = std::make_unique<Solver>();

    Eigen::Vector4d theta_prop = Eigen::Vector4d::Zero();
    // Matrices such that L^(-1) * P * b_i with L cholesky factor of Q_beta
    Eigen::MatrixXd W_beta_curr, W_beta_prop;
=======
    SpMat Q_curr;
    SpMat buffer_1ring, buffer_2ring;

    std::unique_ptr<Solver> Q_solver   = std::make_unique<Solver>();
>>>>>>> main

    // Buffer used to not allocate memory in the steps
    Eigen::VectorXd buffer_V;
    Eigen::MatrixXd buffer_VM;

    std::mt19937_64 rng;
<<<<<<< HEAD
    bool new_theta = false;

    /**
     * @brief Assembles the anisotropic GMRF precision matrix \f$ Q(\theta) \f$ on the triangulated sphere.
     *
     * Evaluates the half-angle diffusion coefficients \f$ c_0(v) = \cosh(\|v\|) \f$ and
     * \f$ c_{1,2}(v) = v_{1,2} \frac{\sinh(\|v\|)}{\|v\|} \f$ (using a second-order Taylor expansion
     * for \f$ \|v\| < 10^{-8} \f$ to avoid removable division-by-zero), builds
     * \f$ K(\theta) = \kappa^2 \tilde{C} + c_0 G_x + c_1 G_y + c_2 G_{xy} \f$ inside `mem.buffer_1ring`,
     * and computes:
     * \f[
     *   Q(\theta) = \frac{1}{4\pi \kappa^2 \sigma_u^2} K(\theta) \tilde{C}^{-1} K(\theta).
     * \f]
     *
     * @param[in]     context Read-only sampler context containing \f$ \tilde{C}, \tilde{C}^{-1}, G_x, G_y, G_{xy} \f$.
     * @param[in,out] mem     Chain workspace whose `buffer_1ring` is overwritten with \f$ K(\theta) \f$.
     * @param[in]     theta   Hyperparameter vector \f$ (\kappa, \sigma_u, v_1, v_2)^\top \f$.
     * @param[out]    out     Destination sparse matrix populated with \f$ Q(\theta) \f$ (2-ring pattern).
     *
     * @note Leaves `mem.buffer_1ring` dirty with \f$ K(\theta) \f$ to avoid redundant zeroing overhead.
     */
    SamplerMemory(const SamplerContext& context, uint_fast64_t seed);
};


=======
    bool new_theta = true;
    size_t curr_iter = 0; 
    SamplerMemory(const SamplerContext& context, std::seed_seq& seed) 
    : buffer_V(context.num_vertices), buffer_VM(context.num_vertices, context.num_years),rng(seed){
        buffer_1ring = context.Gx + context.Gy + context.Gxy;
        buffer_1ring.diagonal() += context.C;
        buffer_2ring = buffer_1ring*context.C_inv.asDiagonal()*buffer_1ring;
        Q_solver->analyzePattern(buffer_2ring);
    }
};

>>>>>>> main
/**
 * @struct ChainOutput
 * @brief Container storing thinned post-burn-in posterior draws for a single chain.
 *
 * Uses a Structure-of-Arrays (SoA) layout mapped by parameter name for direct export to R.
 */
struct ChainOutput{
<<<<<<< HEAD
    std::map<std::string,Eigen::VectorXd> scalars;
    std::map<std::string,Eigen::MatrixXd> vectors;
    std::map<std::string,Eigen::MatrixXd> matrices; // Recorded flattened matrices
=======
    std::map<std::string,std::vector<double>> scalars;
    std::map<std::string,std::vector<Eigen::VectorXd>> vectors;
    std::map<std::string,std::vector<Eigen::MatrixXd>> matrices;
>>>>>>> main
};

/**
 * @brief Assembles the anisotropic GMRF precision matrix \f$ Q(\theta) \f$ on the triangulated sphere.
 *
<<<<<<< HEAD
 * Evaluates the half-angle diffusion coefficients \f$ c_0(v) = \cosh(\|v\|) \f$ and
 * \f$ c_{1,2}(v) = v_{1,2} \frac{\sinh(\|v\|)}{\|v\|} \f$ (using a second-order Taylor expansion
 * for \f$ \|v\| < 10^{-8} \f$ to avoid removable division-by-zero), builds
 * \f$ K(\theta) = \kappa^2 \tilde{C} + c_0 G_x + c_1 G_y + c_2 G_{xy} \f$ inside `mem.buffer_1ring`,
 * and computes:
 * \f[
 *   Q(\theta) = \frac{1}{4\pi \kappa^2 \sigma_u^2} K(\theta) \tilde{C}^{-1} K(\theta).
 * \f]
 *
=======
>>>>>>> main
 * @param[in]     context Read-only sampler context containing \f$ \tilde{C}, \tilde{C}^{-1}, G_x, G_y, G_{xy} \f$.
 * @param[in,out] mem     Chain workspace whose `buffer_1ring` is overwritten with \f$ K(\theta) \f$.
 * @param[in]     theta   Hyperparameter vector \f$ (\kappa, \sigma_u, v_1, v_2)^\top \f$.
 * @param[out]    out     Destination sparse matrix populated with \f$ Q(\theta) \f$ (2-ring pattern).
 *
 * @note Leaves `mem.buffer_1ring` dirty with \f$ K(\theta) \f$ to avoid redundant zeroing overhead.
 */
inline void build_spde_matrix(
    const SamplerContext& context,
    SamplerMemory& mem,
    const Eigen::Vector4d& theta, // (kappa, sigma_u, v1, v2)
    Eigen::SparseMatrix<double>& out
) {
    const double kappa   = theta[ThetaIdx::KAPPA];
    const double sigma_u = theta[ThetaIdx::SIGMA_BETA];
    const double v1      = theta[ThetaIdx::V1];
    const double v2      = theta[ThetaIdx::V2];

    const double v_norm = std::hypot(v1, v2);
    const double c0 = std::cosh(v_norm);
    const double s  = (v_norm < 1e-8)
                    ? (1.0 + (v_norm * v_norm) / 6.0)
                    : (std::sinh(v_norm) / v_norm);

<<<<<<< HEAD
    // 1. K = c0*Gx + c1*Gy + c2*Gxy + kappa^2 * diag(C)
    mem.buffer_1ring = c0 * context.Gx + (v1 * s) * context.Gy + (v2 * s) * context.Gxy;
    mem.buffer_1ring.diagonal() += (kappa * kappa) * context.C;
=======
    const double h11 = c0 + v1 * s;
    const double h22 = c0 - v1 * s;
    const double h12 = v2 * s;

    mem.buffer_1ring = h11 * context.Gx 
    + h22 * context.Gy 
    + h12 * context.Gxy;
    mem.buffer_1ring.diagonal() += (kappa*kappa) * context.C;
>>>>>>> main

    // 2. Q = scale * K * C^{-1} * K
    const double scale = 1.0 / (4.0 * EIGEN_PI * kappa * kappa * sigma_u * sigma_u);
    out = scale * (mem.buffer_1ring * context.C_inv.asDiagonal() * mem.buffer_1ring);
<<<<<<< HEAD
}




inline SamplerMemory::SamplerMemory(const SamplerContext& context, uint_fast64_t seed)
        : W_beta_curr(context.num_vertices, context.num_years),
          W_beta_prop(context.num_vertices, context.num_years),
          buffer_V(context.num_vertices),
          buffer_VM(context.num_vertices, context.num_years),
          rng(seed) {
            const Eigen::Vector4d dummy_theta(1.0,1.0,0.5,0.5);
            build_spde_matrix(context,*this,dummy_theta,buffer_2ring);


            buffer_1ring = context.PhiTPhi;
            buffer_1ring.diagonal().array() += 1.0; // Fallback to have a definite positive matrix

            mu_solver->analyzePattern(buffer_2ring);
            asis_solver->analyzePattern(buffer_1ring);
            mh_curr->analyzePattern(buffer_2ring);
            beta_curr->analyzePattern(buffer_2ring);
            mh_prop->analyzePattern(buffer_2ring);
            beta_prop->analyzePattern(buffer_2ring);
    }
=======
}
>>>>>>> main
