#include "GibbsSteps.hpp"
#include "SamplingUtils.hpp"
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

void Sigma2Step::execute(
    const SamplerContext &ctx,
    ModelState &state,
    SamplerMemory &mem)
{
    // build the sum of squared residuals
    // B * phiTY
    const double cross_term = state.spatial_coeffs.cwiseProduct(ctx.PhiTY).sum();
    // B * PhiTPhi * B
    mem.buffer_VM.noalias() = ctx.PhiTPhi * state.spatial_coeffs; // W = PhiTPhi · B
    const double quad_term = state.spatial_coeffs.cwiseProduct(mem.buffer_VM).sum();
    // sum os squared residuals
    const double sum_sqr_res = std::max(0.0, ctx.YTY_sum - 2.0 * cross_term + quad_term);

    // posterior parameters
    double a_post = a_prior + 0.5 * ctx.num_years * ctx.num_observations;
    double b_post = b_prior + 0.5 * sum_sqr_res;

    // sampling
    std::gamma_distribution<double> gamma(a_post, 1.0 / b_post); // scale = 1/rate
    double sigma2 = 1.0 / gamma(mem.rng);
    state.error_variance = sigma2;
}

void MuCenteredStep::batch_matrix(
    const SamplerContext &ctx,
    const ModelState &state,
    SamplerMemory &mem,
    Factorizer &fact)
{
    if (!mem.new_theta)
        return;

    // Computing the posterior Q_mu matrix
    matrix = ctx.num_years * mem.Q_curr;
    matrix.diagonal().array() += prec_prior;

    fact.submit_matrix(solver, matrix);
}

void MuCenteredStep::execute(
    const SamplerContext &ctx,
    ModelState &state,
    SamplerMemory &mem)
{
    // build the rhs: (1/sigma_mu^2) * mu_prior + Q(theta) * sum_i beta_i
    mem.buffer_V.noalias() = mem.Q_curr * state.spatial_coeffs.rowwise().sum();
    mem.buffer_V.noalias() += prec_prior * mu_prior;

    state.gp_mean = fast_sampling_GMRF(mem.buffer_V, solver, mem.rng);
};

void AsisStep::execute(
    const SamplerContext &ctx,
    ModelState &state,
    SamplerMemory &mem)
{
    mem.buffer_VM = state.spatial_coeffs.colwise() - state.gp_mean;

    mem.buffer_V.noalias() = ctx.PhiTY_sum - ctx.PhiTPhi * mem.buffer_VM.rowwise().sum();
    mem.buffer_V.noalias() += mu_prior * prior_prec_mu;

    state.gp_mean = fast_sampling_GMRF(mem.buffer_V, solver, mem.rng);
    state.spatial_coeffs = mem.buffer_VM.colwise() + state.gp_mean;
};

void AsisStep::batch_matrix(
    const SamplerContext &ctx,
    const ModelState &state,
    SamplerMemory &mem,
    Factorizer &fact)
{
    matrix = (ctx.num_years / state.error_variance) * ctx.PhiTPhi;
    matrix.diagonal().array() += 1.0 / sigma_mu2;

    fact.submit_matrix(solver, matrix);
}

void BlockMHStep::batch_matrix(
    const SamplerContext &ctx,
    const ModelState &state,
    SamplerMemory &mem,
    Factorizer &fact)
{
    adapter->update(state.spde_hyperparameters, last_accept_prob); //! Adapter called in batch_matrix, better to pass hyperparameters not in log form
    build_spde_matrix(ctx, mem, theta_prop, Q_beta_prop);

    // Computing the posterior for our beta_current, and the beta with the proposal
    Q_beta_curr = mem.Q_curr + (1.0 / state.error_variance) * ctx.PhiTPhi;
    Q_beta_prop = Q_prop + (1.0 / state.error_variance) * ctx.PhiTPhi;

    fact.submit_matrix(*Q_solver_prop, Q_prop);
    fact.submit_matrix(*B_solver_curr, Q_beta_curr);
    fact.submit_matrix(*B_solver_prop, Q_beta_prop);
}

void BlockMHStep::execute(
    const SamplerContext &ctx,
    ModelState &state,
    SamplerMemory &mem)
{

    // 1. collapsed log-targets at current and proposed theta, with the current mu, sigma2.
    double log_trgt_curr = eval_log_target(
        state.spde_hyperparameters,
        mem.Q_curr,
        *mem.Q_solver,
        *B_solver_curr,
        W_beta_curr,
        state,
        ctx,
        mem);
    double log_trgt_prop = eval_log_target(
        theta_prop,
        Q_prop, 
        *Q_solver_prop, 
        *B_solver_prop, 
        W_beta_prop, 
        state, 
        ctx, 
        mem);

    // 2. accept/reject.
    double log_a = log_trgt_prop - log_trgt_curr;
    last_accept_prob = std::min(1.0, std::exp(log_a));

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool accept = std::log(unif(mem.rng)) < log_a;

    if (accept)
    {
        mem.new_theta = true;
        state.spde_hyperparameters = theta_prop;
        // Downstream steps and the beta draw below read mem.Q_spde_curr, so move the
        // accepted proposal matrix into the current slot
        mem.Q_curr.swap(Q_prop);
        mem.Q_solver.swap(Q_solver_prop);
        Q_beta_curr.swap(Q_beta_prop);
        B_solver_curr.swap(B_solver_prop);
        W_beta_curr.swap(W_beta_prop);
    }
    else
        mem.new_theta = false;

    // 3. beta draw, conditioned on the accepted theta.
    // W_beta_curr already holds L^{-1} P b_post. Add standard-normal noise:
    std::normal_distribution<double> std_normal(0.0, 1.0);
    for (int j = 0; j < ctx.num_years; ++j)
        for (int i = 0; i < ctx.num_vertices; ++i)
            (W_beta_curr)(i, j) += std_normal(mem.rng);

    // w = L^{-T} (L^{-1} P b_post + z)
    B_solver_curr->matrixU().solveInPlace(W_beta_curr);

    // x = P^T w
    state.spatial_coeffs.noalias() = B_solver_curr->permutationPinv() * W_beta_curr;

    // 4. feed the realized position to the proposal adapter.
    // This refreshes the proposal Cholesky that next sweep's
    // recipe will read to sample theta_prop
    Eigen::Vector4d eta;
    eta(0) = std::log(state.spde_hyperparameters[0]);
    eta(1) = std::log(state.spde_hyperparameters[1]);
    eta(2) = state.spde_hyperparameters[2];
    eta(3) = state.spde_hyperparameters[3];
};

double BlockMHStep::eval_log_target(
    const Eigen::Vector4d &theta,
    const SpMat &Q,
    const Solver &Q_factor,
    const Solver &Q_b_factor,
    Eigen::MatrixXd &W_beta,
    const ModelState &state,
    const SamplerContext &ctx,
    SamplerMemory &mem)
{

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
    W_beta = Q_b_factor.permutationP() * mem.buffer_VM; // W = P b_post
    Q_b_factor.matrixL().solveInPlace(W_beta);          // W = L^{-1} P b_post

    // compute log determinants
    double log_det_Q = 2.0 * Q_factor.matrixL().nestedExpression().diagonal().array().log().sum();
    double log_det_Q_b_post = 2.0 * Q_b_factor.matrixL().nestedExpression().diagonal().array().log().sum();

    double res = std::log(theta(0)) + std::log(theta(1)) + prior->log_pdf(theta) + ctx.num_years * 0.5 * (log_det_Q - log_det_Q_b_post - state.gp_mean.dot(mem.buffer_V)) + 0.5 * W_beta.squaredNorm();

    return res;
};
