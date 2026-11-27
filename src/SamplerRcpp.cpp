// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <random>
#include <memory>
#include <stdexcept>
#include <omp.h>

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
    
    Eigen::initParallel();
    Eigen::setNbThreads(1); 
    
    omp_set_nested(1);
    omp_set_max_active_levels(2);

    using SpMat = Eigen::SparseMatrix<double>;

    
    const long num_obs = Y.rows();
    const long num_years = Y.cols();
    const long num_vertices = Phi.cols();

    SpMat PhiTPhi = Phi.transpose() * Phi;
    PhiTPhi.makeCompressed(); 

    Eigen::MatrixXd PhiTY = Phi.transpose() * Y;
    double YTY_sum = Y.squaredNorm();

    // 3. Parallel configuration 3 chain x 5 cores in gibbs samplers
    ParallelConfig config;
    config.num_chains = 3;
    config.chol_threads = 5;

    SamplerContext ctx(
        PhiTPhi, PhiTY, C, Gx, Gy, Gxy,
        YTY_sum, num_iters, config, num_obs
    );

    // 4. Multiple chains init states
    std::vector<ModelState> init_states;
    init_states.reserve(config.num_chains);

    // Using a generator to compute different states
    std::mt19937 rng(4242); 
    std::uniform_real_distribution<double> runif_kappa(0.5, 3.0); // Da molto range a poco range
    std::uniform_real_distribution<double> runif_sigma(0.5, 2.0);
    std::uniform_real_distribution<double> runif_v(-0.7, 0.7);    // Anisotropia in varie direzioni
    std::uniform_real_distribution<double> runif_eps(0.01, 1.0);  // Varianza errore

    for (size_t c = 0; c < config.num_chains; ++c) {
        ModelState st;
        st.spatial_coeffs = Eigen::MatrixXd::Zero(num_vertices, num_years);
        st.gp_mean = Eigen::VectorXd::Zero(num_vertices);
        
        st.spde_hyperparameters = Eigen::Vector4d(
            runif_kappa(rng), 
            runif_sigma(rng), 
            runif_v(rng), 
            runif_v(rng)
        );
        st.error_variance = runif_eps(rng);
        
        init_states.push_back(st);
    }

    // 5. Symbolic Factorization
    SpMat I_v(num_vertices, num_vertices);
    I_v.setIdentity();

    SpMat matrix_1ring = Gx + Gy + Gxy + I_v;
    SpMat matrix_2ring = matrix_1ring * matrix_1ring; 
    SpMat asis_pattern = PhiTPhi + I_v;
    SpMat mh_pattern = matrix_2ring + PhiTPhi;

    // 6. ChainRunner
    uint64_t base_seed = 42;
    ChainRunner runner(ctx, init_states, base_seed, num_iters, num_burn_in);

    Progress p(num_iters, true);

    runner.configure([&](GibbsEngine& engine) {
        Eigen::VectorXd mu_prior = Eigen::VectorXd::Constant(num_vertices,15.0);
        double prior_var = 100.0;

        engine.add_step(std::make_unique<MuCenteredStep>(mu_prior, prior_var, matrix_2ring));
        engine.add_step(std::make_unique<AsisStep>(mu_prior, prior_var, asis_pattern));
        
        Eigen::Matrix4d empirical_cov = Eigen::Matrix4d::Identity() * 0.01;

        engine.add_step(std::make_unique<BlockMHStep>(
            std::make_unique<EGPrior>(0.53, 1.20, 1.13),
            std::make_unique<ProposalHaario>(num_burn_in, 100, empirical_cov),
            mh_pattern
        ));
        
        engine.add_step(std::make_unique<Sigma2Step>(1.0, 0.1));

        // Recorder
        engine.add_collector(std::make_unique<ScalarRecorder>(
            num_iters, 0, 1, "new_theta",
            [](const SamplerContext& ctx, const ModelState& s, const SamplerMemory& mem){
                return mem.new_theta ? 1.0 : 0.0;
            }
        ));
        engine.add_collector(std::make_unique<ScalarRecorder>(
            num_iters, num_burn_in, 1, "sigma2_epsilon",
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
            num_iters, num_burn_in, 100, num_vertices, "mu",
            [](const SamplerContext& c, const ModelState& s, const SamplerMemory& m) {
                return s.gp_mean;
            }
        ));

        engine.add_collector(std::make_unique<MatrixRecorder>(
            num_iters,num_burn_in,150,num_vertices,num_years,"beta",
            [](const SamplerContext& c, const ModelState& s, const SamplerMemory& m){
                return s.spatial_coeffs;
            }
        ));

        // Modifier protetto per OpenMP: Rcpp NON può essere chiamato dai background thread
        engine.add_modifier([&p](const SamplerContext& c, const SamplerMemory& m, const ModelState& s) {
            if (omp_get_thread_num() == 0) {
                if (m.curr_iter % 50 == 0) {
                    Rcpp::checkUserInterrupt();
                }
                p.increment();
            }
        });
    });

    // 7. Esecuzione Parallela
    std::vector<ChainOutput> outs = runner.run();

    // 8. Ritorna tutte le catene al backend R
    Rcpp::List r_chains(config.num_chains);
    for (size_t i = 0; i < config.num_chains; ++i) {
        r_chains[i] = Rcpp::List::create(
            Rcpp::Named("sigma_epsilon") = outs[i].scalars.at("sigma2_epsilon"),
            Rcpp::Named("theta")         = outs[i].vectors.at("theta"),
            Rcpp::Named("mu")            = outs[i].vectors.at("mu"),
            Rcpp::Named("beta")          = outs[i].matrices.at("beta"),
            Rcpp::Named("theta_accept")  = outs[i].scalars.at("new_theta")
        );
    }

    return r_chains;
}