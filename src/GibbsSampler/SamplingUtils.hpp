/**
 - @file sampler_types.hpp
 - @brief Core data structures, preallocated MCMC workspace, and SPDE precision matrix assembly.
 *
 - Defines the state, read-only sufficient statistics, thread-local solver memory, and output containers.
 */
#pragma once

#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>


/**
 - @enum ThetaIdx
 - @brief Positional indices for the SPDE hyperparameter vector \f$ \theta = (\kappa, \sigma_u, v_1, v_2) \f$.
 */
enum ThetaIdx : int {
  KAPPA      = 0,
    SIGMA_BETA = 1,
    V1         = 2,
    V2         = 3
};

/**
 - @struct ParallelConfig
 - @brief Threading configuration across MCMC chains and linear algebra routines.
 */
struct ParallelConfig {
  int num_chains = 1;
  int chol_threads = 1;
};


/**
 - @struct ModelState
 - @brief Mutable state of a single Markov chain at the current Gibbs iteration.
 */
struct ModelState {
  Eigen::MatrixXd spatial_coeffs;
  Eigen::VectorXd spde_hyperparameters;
  Eigen::VectorXd gp_mean;
  double error_variance;
};


/**
 - @struct SamplerContext
 - @brief Immutable, read-only context shared safely across all parallel MCMC chains.
 *
 - Stores precomputed sufficient statistics (eliminating \f$ \mathcal{O}(N) \f$ grid operations inside
 - the MCMC loop), finite element geometric matrices on \f$ \mathbb{S}^2 \f$, and problem dimensions.
 */
struct SamplerContext {
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
  Identity_V([this]() {
    Eigen::SparseMatrix<double> id(num_vertices, num_vertices);
    id.setIdentity();
    return id;
  }()),
  num_iters(num_iters_),
  config(config_) {}

  const Eigen::SparseMatrix<double> PhiTPhi;
  const Eigen::MatrixXd PhiTY;
  const Eigen::VectorXd PhiTY_sum; // rowwise sum of PhiTY

  const Eigen::VectorXd C;     // Since C is diagonal we store it as a vector
  const Eigen::VectorXd C_inv; // We also store the inverse of C

  const Eigen::SparseMatrix<double> Gx;
  const Eigen::SparseMatrix<double> Gy;
  const Eigen::SparseMatrix<double> Gxy;

  const double YTY_sum;

  const Eigen::Index num_observations;
  const Eigen::Index num_vertices;
  const Eigen::Index num_years;

  const Eigen::SparseMatrix<double> Identity_V;

  const size_t num_iters;

  const ParallelConfig config;
};

/**
 - @struct SamplerMemory
 - @brief Preallocated, thread-local workspace for a single MCMC chain.
 *
 - Eliminates dynamic heap allocations during Gibbs sweeps by holding reusable sparse
 - buffers, dense scratchpads, and symbolically analyzed sparse Cholesky solvers grouped
 - into 2-ring (SPDE / posterior latent field) and 1-ring (non-centered ASIS) sparsity patterns.
 */
struct SamplerMemory {
  using SpMat  = Eigen::SparseMatrix<double>;
  using Solver = Eigen::SimplicialLLT<SpMat>;

  SpMat Q_curr;
  SpMat buffer_1ring, buffer_2ring;

  std::unique_ptr<Solver> Q_solver = std::make_unique<Solver>();

  // Buffer used to not allocate memory in the steps
  Eigen::VectorXd buffer_V;
  Eigen::MatrixXd buffer_VM;

  std::mt19937_64 rng;
  bool new_theta = true;
  size_t curr_iter = 0;

  SamplerMemory(const SamplerContext& context, std::seed_seq& seed)
    : buffer_V(context.num_vertices),
      buffer_VM(context.num_vertices, context.num_years),
      rng(seed) {
    buffer_1ring = context.Gx + context.Gy + context.Gxy;
    buffer_1ring.diagonal() += context.C;
    buffer_2ring = buffer_1ring * context.C_inv.asDiagonal() * buffer_1ring;
    Q_solver->analyzePattern(buffer_2ring);
  }
};

/**
 - @struct ChainOutput
 - @brief Container storing thinned post-burn-in posterior draws for a single chain.
 *
 - Uses a Structure-of-Arrays (SoA) layout mapped by parameter name for direct export to R.
 */
struct ChainOutput {
  std::map<std::string, std::vector<double>> scalars;
  std::map<std::string, std::vector<Eigen::VectorXd>> vectors;
  std::map<std::string, std::vector<Eigen::MatrixXd>> matrices;
};

/**
 - @brief Assembles the anisotropic GMRF precision matrix \f$ Q(\theta) \f$ on the triangulated sphere.
 *
 - @param[in]     context Read-only sampler context containing \f$ \tilde{C}, \tilde{C}^{-1}, G_x, G_y, G_{xy} \f$.
 - @param[in,out] mem     Chain workspace whose `buffer_1ring` is overwritten with \f$ K(\theta) \f$.
 - @param[in]     theta   Hyperparameter vector \f$ (\kappa, \sigma_u, v_1, v_2)^\top \f$.
 - @param[out]    out     Destination sparse matrix populated with \f$ Q(\theta) \f$ (2-ring pattern).
 *
 - @note Leaves `mem.buffer_1ring` dirty with \f$ K(\theta) \f$ to avoid redundant zeroing overhead.
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

  const double h11 = c0 + v1 * s;
  const double h22 = c0 - v1 * s;
  const double h12 = v2 * s;

  mem.buffer_1ring = h11 * context.Gx
  + h22 * context.Gy
  + h12 * context.Gxy;
  mem.buffer_1ring.diagonal() += (kappa * kappa) * context.C;

  // 2. Q = scale * K * C^{-1} * K
  const double scale = 1.0 / (4.0 * EIGEN_PI * kappa * kappa * sigma_u * sigma_u);
  out = scale * (mem.buffer_1ring * context.C_inv.asDiagonal() * mem.buffer_1ring);
}
