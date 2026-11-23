#include "GibbsSteps.hpp"
#include "SamplingUtils.hpp"
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

void Sigma2Step::execute(ModelState &state,
                         const SamplerContext &ctx,
                         SamplerMemory &mem)
{
    // build the sum of squared residuals
    // B * phiTY
    const double cross = state.spatial_coeffs.cwiseProduct(ctx.PhiTY).sum();
    // B * PhiTPhi * B
    mem.buffer_VM.noalias() = ctx.PhiTPhi * state.spatial_coeffs; // W = PhiTPhi · B
    const double quad = state.spatial_coeffs.cwiseProduct(mem.buffer_VM).sum();
    // sum os squared residuals
    const double ssr = std::max(0.0, ctx.YTY - 2.0 * cross + quad);

    // posterior parameters
    double a_post = a0 + 0.5 * ctx.num_years * ctx.num_observations;
    double b_post = b0 + 0.5 * ssr;

    // sampling
    std::gamma_distribution<double> gamma(a_post, 1.0 / b_post); // scale = 1/rate
    double sigma2 = 1.0 / gamma(mem.rng);
    state.error_variance = sigma2;
}

Eigen::VectorXd fast_sampling_GMRF(
    Eigen::Ref<const Eigen::VectorXd> b,
    const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &Q_factor,
    std::mt19937_64 &rng)
{
    Eigen::VectorXd sample = Q_factor.permutationP() * b;
    sample = Q_factor.matrixL().solve(sample);
    std::normal_distribution<double> std_normal(0.0, 1.0);
    for (int i = 0; i < sample.size(); ++i)
    {
        sample(i) += std_normal(rng);
    }
    Eigen::VectorXd buffer = Q_factor.matrixU().solve(sample);
    return Q_factor.permutationPinv() * buffer;
};

void MuCenteredStep::execute(ModelState &state,
                             const SamplerContext &ctx,
                             SamplerMemory &mem)
{
    // build the rhs: (1/sigma_mu^2) * mu_prior + Q(theta) * sum_i beta_i   (Eq. 24)
    mem.buffer_V.noalias() = mem.Q_spde_curr * state.spatial_coeffs.rowwise().sum();
    mem.buffer_V.noalias() += prior_prec_mu * mu_prior;

    state.gp_mean = fast_sampling_GMRF(mem.buffer_V, solver_, mem.rng);
};

void AsisStep::execute(ModelState &state,
                       const SamplerContext &ctx,
                       SamplerMemory &mem)
{
    mem.buffer_VM = state.spatial_coeffs.colwise() - state.gp_mean;

    mem.buffer_V.noalias() = ctx.PhiTY_sum - ctx.PhiTPhi * mem.buffer_VM.rowwise().sum();
    mem.buffer_V.noalias() += mu_prior * prior_prec_mu;

    state.gp_mean = fast_sampling_GMRF(mem.buffer_V, solver_, mem.rng);
    state.spatial_coeffs = mem.buffer_VM.colwise() + state.gp_mean;
};

double PCPrior::f(double v_norm)
{
    return std::sqrt(EIGEN_PI * (3 * std::cosh(2 * v_norm) + 1) / 3);
};

double PCPrior::f_prime(double v_norm)
{
    return std::sqrt(EIGEN_PI) * std::sinh(2 * v_norm) + 1 / std::sqrt(std::cosh(2 * v_norm) + 1 / 3);
};

double PCPrior::f_prime_over_v(double v_norm)
{
    // Taylor expansion of third order of f'(v)/v
    // centered in 0 is f''(0)
    // f''(0) = lim_{v->0} f'(v)/v = sqrt(3*pi)
    static const double f2_0 = std::sqrt(3.0 * EIGEN_PI);
    if (v_norm < 1e-8)
    {
        return f2_0;
    }
    else
    {
        return f_prime(v_norm) / v_norm;
    }
};

double PCPrior::log_pdf(const Eigen::Vector4d &theta) const
{
    double v_norm = std::sqrt(theta(2) * theta(2) + theta(3) * theta(3));

    double res = log(lambda_kappa) + log(lambda_v) + log(f(v_norm)) + log(f_prime_over_v(v_norm)) - log(2 * EIGEN_PI) - lambda_v * (f(v_norm) - f(0)) - lambda_kappa * f(v_norm) * theta(0) + log(lambda_sigma) - theta(1) * lambda_sigma;

    return res;
};

double MHStepBlockThetaBeta::eval_collapsed_log_target(const Eigen::Vector4d &theta,
                                                       const ModelState &state,
                                                       const SamplerContext &ctx,
                                                       SamplerMemory &mem,
                                                       bool use_prop)
{
    // choose the pre-assembled matrix and pre-computed factors for this slot
    SpMat &Q = use_prop ? mem.Q_spde_prop : mem.Q_spde_curr;
    Solver &Q_factor = use_prop ? Q_solver_prop : Q_solver_curr;
    Solver &B_factor = use_prop ? B_solver_prop : B_solver_curr;

    Eigen::MatrixXd &W = W_beta_scratch;

    // Q * mu
    mem.buffer_V = Q * state.gp_mean;

    // b_i,post for all years at once (V x M)
    // b_i,post = Q(theta) * mu + Phi^T Y_i / sigma^2
    mem.buffer_VM = ctx.PhiTY / state.error_variance;
    mem.buffer_VM.colwise() += mem.buffer_V;

    // Compute sum_i b_i^T Q_beta,post^{-1} b_i.
    // Eigen factors: P Q_beta,post P^T = L L^T (P^{-1}=P^T),
    // so Q_beta,post^{-1} = P^T L^{-T} L^{-1} P  =>  b_i^T Q_beta,post^{-1} b_i = ||L^{-1} P b_i||^2.
    // Stack all years into W and the whole sum is W.squaredNorm() (Frobenius)
    // W needs its own buffer since the permutation writes a fresh matrix
    W = B_factor.permutationP() * mem.buffer_VM; // W = P b_post
    B_factor.matrixL().solveInPlace(W);          // W = L^{-1} P b_post

    // compute log determinants
    double log_det_Q = 2.0 * Q_factor.matrixL().nestedExpression().diagonal().array().log().sum();
    double log_det_Q_b_post = 2.0 * B_factor.matrixL().nestedExpression().diagonal().array().log().sum();

    double res = std::log(theta(0)) + std::log(theta(1)) + prior->log_pdf(theta) + ctx.num_years * 0.5 * (log_det_Q - log_det_Q_b_post - state.gp_mean.dot(mem.buffer_V)) + 0.5 * W.squaredNorm();

    return res;
};

void MHStepBlockThetaBeta::execute(ModelState &state,
                                   const SamplerContext &ctx,
                                   SamplerMemory &mem)
{

    // 1. collapsed log-targets at current and proposed theta, with the current mu, sigma2.
    double log_trgt_curr = eval_collapsed_log_target(state.spde_hyperparameters,
                                                     state, ctx, mem, 0);
    double log_trgt_prop = eval_collapsed_log_target(mem.theta_prop,
                                                     state, ctx, mem, 1);

    // 2. accept/reject.
    double log_a = log_trgt_prop - log_trgt_curr;
    last_accept_prob = std::min(1.0, std::exp(log_a));

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool accept = std::log(unif(mem.rng)) < log_a;

    if (accept)
    {
        mem.new_theta = true;
        state.spde_hyperparameters = mem.theta_prop;
        // Downstream steps and the beta draw below read mem.Q_spde_curr, so move the
        // accepted proposal matrix into the current slot
        // The factors are not swapped since they are rebuilt next sweep by the recipe
        mem.Q_spde_curr.swap(mem.Q_spde_prop);
    }
    else
    {
        mem.new_theta = false;
    }

    // 3. beta draw, conditioned on the accepted theta.
    // The accepted Q_beta,post factor is the proposal one iff we accepted
    Solver &B_acc = accept ? B_solver_prop : B_solver_curr;

    // rhsides: PhiTY / sigma2 + Q(theta_acc) * mu
    mem.buffer_VM = ctx.PhiTY / state.error_variance;
    mem.buffer_VM.colwise() += mem.Q_spde_curr * state.gp_mean;

    // draw one column (year) at a time
    for (int i = 0; i < ctx.num_years; ++i)
    {
        state.spatial_coeffs.col(i) = fast_sampling_GMRF(mem.buffer_VM.col(i), B_acc, mem.rng);
    }

    // 4. feed the realized position to the proposal adapter.
    // This refreshes the proposal Cholesky that next sweep's
    // recipe will read to sample theta_prop
    Eigen::Vector4d eta;
    eta(0) = std::log(state.spde_hyperparameters[0]);
    eta(1) = std::log(state.spde_hyperparameters[1]);
    eta(2) = state.spde_hyperparameters[2];
    eta(3) = state.spde_hyperparameters[3];
    adapter->update(eta, last_accept_prob);
};
