#pragma once

#include <Eigen/Dense>
#include "SamplingUtils.hpp"
#include "StepUtils.hpp"



class Factorizer{
using SpMat = Eigen::SparseMatrix<double>;
using Solver = Eigen::SimplicialLLT<SpMat>;
public:
    Factorizer(size_t num_threads_) : num_threads(num_threads_){
        batch.reserve(10);
    }

    void submit_matrix(Solver& solver,const SpMat& matrix){
        batch.push_back({&solver,&matrix});
    }

    void factorize_parallel(){
        //TODO Add openMP directive
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
        std::unique_ptr<ProposalAdapter> adapter_,
        const SpMat& matrix_2ring
    ) : prior(std::move(prior_)), 
    adapter(std::move(adapter_)),
    Q_prop(matrix_2ring),
    Q_beta_curr(matrix_2ring),
    Q_beta_prop(matrix_2ring){
        Q_solver_prop.analyzePattern(Q_prop);
        B_solver_prop.analyzePattern(Q_prop);
        B_solver_curr.analyzePattern(Q_prop);
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
    double eval_collapsed_log_target(
        const Eigen::Vector4d &theta,
        const ModelState &state,
        const SamplerContext &ctx,
        SamplerMemory &mem,
        bool use_prop
    );


    std::unique_ptr<ThetaPrior> prior;
    std::unique_ptr<ProposalAdapter> adapter;
    double last_accept_prob = 0.0; // consumed only by the adapter


    SpMat Q_prop, Q_beta_prop,Q_beta_curr;

    // Factor of Q(theta) for the proposal slot
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


    Eigen::Vector4d theta_prop;

};
