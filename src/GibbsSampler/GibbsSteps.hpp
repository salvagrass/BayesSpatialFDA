#pragma once

#include <Eigen/Dense>
#include "SamplingUtils.hpp"
#include "StepUtils.hpp"
#include <omp.h>


class Factorizer{
using SpMat = Eigen::SparseMatrix<double>;
using Solver = Eigen::SimplicialLLT<SpMat>;
public:
    Factorizer(const ParallelConfig& config) : num_threads(config.chol_threads){
        batch.reserve(10);
    }

    void submit_matrix(Solver& solver,const SpMat& matrix){
        batch.push_back({&solver,&matrix});
    }

    void factorize_parallel(){

        #pragma omp parallel for num_threads(num_threads) if (num_threads > 1) schedule(static)
        for(size_t i = 0; i<batch.size();i++){
            batch[i].solver->factorize(*batch[i].matrix);
        }

        batch.clear();
    }

private:
    struct batch_element{
        Solver* solver;
        const SpMat* matrix;
    };
    std::vector<batch_element> batch;
    size_t num_threads;
};


class Step{
public:
    Step() = default;
    virtual ~Step() = default;

    virtual void execute(
        const SamplerContext& ctx,
        ModelState& state,
        SamplerMemory& memory
    ) = 0;

    // A step may contribute a factorization to the sweep-start barrier
    // Steps with no factorization (Sigma2Step) inherit this no-op
    virtual void batch_matrix(
        const SamplerContext& ctx,
        const ModelState& state,
        SamplerMemory& mem,
        Factorizer& factorizer
    ) {return;}

protected:
    using SpMat = Eigen::SparseMatrix<double>;
    using Solver = Eigen::SimplicialLLT<SpMat>;
    //TODO Maybe we could wrap the package in the same namespace
};

class Sigma2Step : public Step{
public:
    Sigma2Step(double a0, double b0) : Step(), a_prior(a0), b_prior(b0) {}
    void execute(
        const SamplerContext& ctx,
        ModelState& state,
        SamplerMemory& mem
    ) override;

private:
    double a_prior;
    double b_prior;

};

class MuCenteredStep : public Step{
public:
    using Step::Step;
    MuCenteredStep(Eigen::VectorXd mu_prior_, double var_prior_,const SpMat& matrix_2ring) 
    : mu_prior(mu_prior_), 
    var_prior(var_prior_), 
    prec_prior(1.0 / var_prior),
    matrix(matrix_2ring){
        solver.analyzePattern(matrix);
    }
    void execute(
        const SamplerContext& ctx,
        ModelState& state,
        SamplerMemory& mem
    ) override;
    void batch_matrix(
        const SamplerContext& ctx,
        const ModelState& state,
        SamplerMemory& mem,
        Factorizer& factorizer
    ) override;

private:
    Eigen::VectorXd mu_prior;
    double var_prior;
    double prec_prior;

    SpMat matrix;
    Solver solver;
};

class AsisStep : public Step{
public:
    using Step::Step;
    AsisStep(Eigen::VectorXd mu_prior_, double sigma_mu2_,const SpMat& matrix_1ring) 
    : mu_prior(mu_prior_), 
    sigma_mu2(sigma_mu2_), 
    prior_prec_mu(1.0 / sigma_mu2_), 
    matrix(matrix_1ring){
        solver.analyzePattern(matrix);
    }
    void execute(
        const SamplerContext& ctx,
        ModelState& state,
        SamplerMemory& mem
    ) override;
    void batch_matrix(
        const SamplerContext& ctx,
        const ModelState& state,
        SamplerMemory& mem,
        Factorizer& factorizer
    ) override;
private:
    Eigen::VectorXd mu_prior;
    double sigma_mu2;
    double prior_prec_mu;

    SpMat matrix;
    Solver solver;
};



class BlockMHStep : public Step{
public:
    BlockMHStep(
        std::unique_ptr<ThetaPrior> prior_,
        std::unique_ptr<ProposalStrategy> proposal_,
        const SpMat& matrix_2ring
    ) : prior(std::move(prior_)), 
    proposal(std::move(proposal_)),
    Q_prop(matrix_2ring),
    Q_beta_prop(matrix_2ring),
    Q_beta_curr(matrix_2ring){
        Q_solver_prop -> analyzePattern(Q_prop);
        B_solver_prop -> analyzePattern(Q_prop);
        B_solver_curr -> analyzePattern(Q_prop);
    }

    void execute(
        const SamplerContext& ctx,
        ModelState& state,
        SamplerMemory& mem
    ) override;
    
    void batch_matrix(
        const SamplerContext& ctx,
        const ModelState& state,
        SamplerMemory& mem,
        Factorizer& factorizer
    ) override;


private:
    // Evaluates log pi(eta) at (theta; current mu, current sigma2 from state)
    double eval_log_target(
        const Eigen::Vector4d &theta,
         const SpMat& Q,
        const Solver& Q_factor, 
        const Solver& Q_b_factor,
        Eigen::MatrixXd& W_beta,
        const ModelState &state,
        const SamplerContext &ctx,
        SamplerMemory &mem
    );

    Eigen::Vector4d theta_to_eta(const Eigen::Vector4d& theta);
    Eigen::Vector4d eta_to_theta(const Eigen::Vector4d& eta);

    std::unique_ptr<ThetaPrior> prior;
    std::unique_ptr<ProposalStrategy> proposal;

    SpMat Q_prop, Q_beta_prop,Q_beta_curr;

    // Factor of Q(theta) for the proposal slot
    std::unique_ptr<Solver> Q_solver_prop = std::make_unique<Solver>();
    // Factors of Q_beta,post(theta) for the current / proposal slots
    //Solver B_solver_curr;
    std::unique_ptr<Solver> B_solver_curr = std::make_unique<Solver>();
    std::unique_ptr<Solver> B_solver_prop = std::make_unique<Solver>();

    // Private assembly scratch for the factor sum_i b_i^T Q_beta,post^{-1} b_i,
    Eigen::MatrixXd W_beta_curr;
    Eigen::MatrixXd W_beta_prop;

    Eigen::Vector4d theta_prop;

};
