#pragma once

#include <Eigen/Dense>
#include "SamplingUtils.hpp"
#include "GibbsHooks.hpp"

using SpMat = Eigen::SparseMatrix<double>;
using Solver = Eigen::SimplicialLLT<SpMat>;

class Factorizer;

// Fast sampling of Gaussian Markov random fields, Rue 2001
Eigen::VectorXd fast_sampling_GMRF(Eigen::Ref<const Eigen::VectorXd> b,
                                   const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &Q_factor,
                                   std::mt19937_64 &rng);

class Step
{
private:
public:
    Step() = default;
    virtual ~Step() = default;
    virtual void execute(ModelState &ModelState_,
                         const SamplerContext &SamplerContext_,
                         SamplerMemory &Samplermemory_) = 0;
    // A step may contribute a factorization to the sweep-start barrier
    // Steps with no factorization (Sigma2Step) inherit this no-op
    virtual void submit_factorization(Factorizer &, ModelState &,
                                      const SamplerContext &,
                                      SamplerMemory &) {}
};

class Sigma2Step : public Step
{
private:
    double a0;
    double b0;

public:
    using Step::Step;
    Sigma2Step(double a0_, double b0_) : Step(), a0(a0_), b0(b0_) {}
    void execute(ModelState &ModelState_,
                 const SamplerContext &SamplerContext_,
                 SamplerMemory &Samplermemory_) override;
};

class MuCenteredStep : public Step
{
private:
    Eigen::VectorXd mu_prior;
    double sigma_mu2;
    double prior_prec_mu;

    SpMat matrix_;
    Solver solver_;

public:
    using Step::Step;
    MuCenteredStep(Eigen::VectorXd mu_prior_, double sigma_mu2_) : Step(), mu_prior(mu_prior_), sigma_mu2(sigma_mu2_), prior_prec_mu(1.0 / sigma_mu2_) {}
    void execute(ModelState &ModelState_,
                 const SamplerContext &SamplerContext_,
                 SamplerMemory &Samplermemory_) override;
    void submit_factorization(Factorizer &f, ModelState &ModelState_,
                              const SamplerContext &ctx, SamplerMemory &mem) override;
};

class AsisStep : public Step
{
private:
    Eigen::VectorXd mu_prior;
    double sigma_mu2;
    double prior_prec_mu;

    SpMat matrix_;
    Solver solver_;

public:
    using Step::Step;
    AsisStep(Eigen::VectorXd mu_prior_, double sigma_mu2_) : Step(), mu_prior(mu_prior_), sigma_mu2(sigma_mu2_), prior_prec_mu(1.0 / sigma_mu2_) {}
    void execute(ModelState &ModelState_,
                 const SamplerContext &SamplerContext_,
                 SamplerMemory &Samplermemory_) override;
    void submit_factorization(Factorizer &f, ModelState &ModelState_,
                              const SamplerContext &ctx,
                              SamplerMemory &mem) override;
};

class ThetaPrior
{
public:
    ThetaPrior() = default;
    virtual ~ThetaPrior() = default;
    virtual double log_pdf(const Eigen::Vector4d &theta) const = 0;
};

class PCPrior : public ThetaPrior
{
private:
    double lambda_kappa, lambda_v, lambda_sigma;
    static double f(double v_norm);
    static double f_prime(double v_norm);
    static double f_prime_over_v(double v_norm);

public:
    using ThetaPrior::ThetaPrior;
    PCPrior(double lk, double lv, double ls) : ThetaPrior(), lambda_kappa(lk), lambda_v(lv), lambda_sigma(ls) {}
    double log_pdf(const Eigen::Vector4d &theta) const override;
};

// Proposal adaptation for the theta block (eta-space, d = 4).
// Supplies the proposal Cholesky L (read by the factorization recipe to form
// eta_prop = eta_curr + L * z) and adapts it from the realized chain.
class ProposalAdapter
{
public:
    ProposalAdapter() = default;
    virtual ~ProposalAdapter() = default;
    // Lower-triangular Cholesky L of the proposal covariance.
    virtual const Eigen::Matrix4d &proposal_chol() const = 0;
    // Adapt after one MH iteration, given the chain's current position (eta)
    // and the acceptance probability just computed.
    virtual void update(const Eigen::Vector4d &eta, double accept_prob) = 0;
};

// Non-adaptive: constant proposal covariance.
class FixedProposal : public ProposalAdapter
{
private:
    Eigen::Matrix4d L_;

public:
    explicit FixedProposal(const Eigen::Matrix4d &L) : L_(L) {}
    const Eigen::Matrix4d &proposal_chol() const override { return L_; }
    void update(const Eigen::Vector4d &, double) override {}
};

class MHStepBlockThetaBeta : public Step
{
private:
    std::unique_ptr<ThetaPrior> prior;
    std::unique_ptr<ProposalAdapter> adapter;
    double last_accept_prob = 0.0; // consumed only by the adapter

    // Factors of Q(theta) for the current / proposal slots
    Solver Q_solver_curr;
    Solver Q_solver_prop;
    // Factors of Q_beta,post(theta) for the current / proposal slots
    Solver B_solver_curr;
    Solver B_solver_prop;
    // Plain values for both since the recipe rebuilds both factorizations at every sweep

    // Private assembly scratch for Q_beta,post before factorizing
    // Used by both curr and prop
    // The B matrix is never read after factorization, so one buffer suffices
    SpMat B_mat_scratch;

    // Private assembly scratch for the factor sum_i b_i^T Q_beta,post^{-1} b_i,
    // One buffer suffices since it is only used inside eval_collapsed_log_target,
    // and those two calls (curr, prop) are sequential in execute
    Eigen::MatrixXd W_beta_scratch;

    // Evaluates log pi(eta) at (theta; current mu, current sigma2 from state)
    double eval_collapsed_log_target(const Eigen::Vector4d &theta,
                                     const ModelState &state,
                                     const SamplerContext &ctx,
                                     SamplerMemory &mem,
                                     bool use_prop);

public:
    using Step::Step;
    MHStepBlockThetaBeta(std::unique_ptr<ThetaPrior> prior_,
                         std::unique_ptr<ProposalAdapter> adapter_)
        : Step(), prior(std::move(prior_)), adapter(std::move(adapter_)) {}
    void execute(ModelState &ModelState_,
                 const SamplerContext &SamplerContext_,
                 SamplerMemory &Samplermemory_) override;
    void submit_factorization(Factorizer &f, ModelState &ModelState_,
                              const SamplerContext &ctx, SamplerMemory &mem) override;
};
