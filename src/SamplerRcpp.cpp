// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp> // <-- Decommentato per la progress bar
#include <random>
#include <memory>
#include <stdexcept>

#include "GibbsSampler/GibbsCore.hpp"
#include "GibbsSampler/SamplingUtils.hpp"
#include "GibbsSampler/GibbsHooks.hpp"
#include "GibbsSampler/GibbsSteps.hpp"
#include "GibbsSampler/StepUtils.hpp"

// [[Rcpp::export(name = ".run_gibbs_sampler_cpp")]]
Rcpp::List run_gibbs_sampler(
    const Eigen::SparseMatrix<double>& Phi,
    const Eigen::MatrixXd& Y,
    const Eigen::SparseMatrix<double>& Gx,
    const Eigen::SparseMatrix<double>& Gy,
    const Eigen::SparseMatrix<double>& Gxy,
    const Eigen::VectorXd& C,
    size_t num_iters,
    size_t num_burn_in
){
    using SpMat = Eigen::SparseMatrix<double>;

    // 1. Statistiche Sufficienti e Dimensioni
    const long num_obs = Y.rows();
    const long num_years = Y.cols();
    const long num_vertices = Phi.cols();

    SpMat PhiTPhi = Phi.transpose() * Phi;
    PhiTPhi.makeCompressed(); 

    Eigen::MatrixXd PhiTY = Phi.transpose() * Y;
    double YTY_sum = Y.squaredNorm();

    // 2. Configurazione Parallela (Singola Catena)
    ParallelConfig config;
    config.num_chains = 1;
    config.chol_threads = 1;
    config.year_threads = 1;

    // 3. Costruzione SamplerContext
    SamplerContext ctx(
        PhiTPhi, PhiTY, C, Gx, Gy, Gxy,
        YTY_sum, num_iters, config, num_obs
    );

    // 4. Inizializzazione ModelState
    ModelState init_state;
    init_state.spatial_coeffs = Eigen::MatrixXd::Zero(num_vertices, num_years);
    init_state.spde_hyperparameters = Eigen::Vector4d(1.0, 1.0, 0.0, 0.0); // kappa, sigma_u, v1, v2
    init_state.gp_mean = Eigen::VectorXd::Zero(num_vertices);
    init_state.error_variance = 1.0;

    std::vector<ModelState> init_states = { init_state };

    // 5. Fattorizzazione Simbolica (Creazione dei Pattern)
    SpMat I_v(num_vertices, num_vertices);
    I_v.setIdentity();

    // Pattern 1-ring e 2-ring per la precisione Q dello SPDE
    SpMat matrix_1ring = Gx + Gy + Gxy + I_v;
    SpMat matrix_2ring = matrix_1ring * matrix_1ring; 

    // Pattern per AsisStep
    SpMat asis_pattern = PhiTPhi + I_v;

    // Pattern per BlockMHStep (Unione pattern spaziale e d'osservazione)
    SpMat mh_pattern = matrix_2ring + PhiTPhi;

    // 6. ChainRunner e Configurazione degli Hook
    uint64_t base_seed = 42;
    ChainRunner runner(ctx, init_states, base_seed, num_iters, num_burn_in);

    // Inizializza la Progress Bar
    Progress p(num_iters, true);

    runner.configure_hook([&](GibbsEngine& engine) {
        Eigen::VectorXd mu_prior = Eigen::VectorXd::Zero(num_vertices);
        double prior_var = 100.0;

        // A. Aggiunta step
        engine.add_step(std::make_unique<MuCenteredStep>(mu_prior, prior_var, matrix_2ring));
        engine.add_step(std::make_unique<AsisStep>(mu_prior, prior_var, asis_pattern));
        engine.add_step(std::make_unique<BlockMHStep>(
            std::make_unique<EGPrior>(1.0, 1.0, 1.0),
            std::make_unique<ProposalHaario>(num_burn_in,100,0.001),
            mh_pattern
        ));
        engine.add_step(std::make_unique<Sigma2Step>(1.0,0.1));

        // B. Aggiunta Recorder

        engine.add_collector(std::make_unique<ScalarRecorder>(
            num_iters,num_burn_in,1,"new_theta",
            [](const SamplerContext& ctx,const ModelState& s,const SamplerMemory& mem){
                return mem.new_theta ? 1 : 0;
            }
        ));
        engine.add_collector(std::make_unique<ScalarRecorder>(
            num_iters,num_burn_in,1,"sigma2_epsilon",
            [](const SamplerContext& ctx, const ModelState& s, const SamplerMemory& mem){
                return s.error_variance;
            }
        ));
        engine.add_collector(std::make_unique<VectorRecorder>(
            num_iters, num_burn_in, 1, 4, "theta",
            [](const SamplerContext& c, const ModelState& s, const SamplerMemory& m) {
                return s.spde_hyperparameters;
            }
        ));
        engine.add_collector(std::make_unique<VectorRecorder>(
            num_iters, num_burn_in, 1, num_vertices, "mu",
            [](const SamplerContext& c, const ModelState& s, const SamplerMemory& m) {
                return s.gp_mean;
            }
        ));

        // C. Modifier per interrompere la catena e aggiornare la progress bar
        engine.add_modifier([&p](const SamplerContext& c, const SamplerMemory& m, const ModelState& s) {
            if (m.curr_iter % 50 == 0) {
                Rcpp::checkUserInterrupt();
            }
            p.increment();
        });
    });

    // 7. Esecuzione
    std::vector<ChainOutput> outs = runner.run();

    return Rcpp::List::create(
        Rcpp::Named("sigma_epsilon") = outs[0].scalars.at("sigma2_epsilon"),
        Rcpp::Named("theta") = outs[0].vectors.at("theta"),
        Rcpp::Named("mu")    = outs[0].vectors.at("mu"),
        Rcpp::Named("theta_accept") = outs[0].scalars.at("new_theta")
    );
}