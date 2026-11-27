#pragma once

#include "SamplingUtils.hpp"
#include "GibbsSteps.hpp"
#include "GibbsHooks.hpp"

#include <omp.h>

class GibbsEngine
{
public:
    GibbsEngine(
        const SamplerContext& ctx, const ModelState& init, 
        std::seed_seq &seed, int num_iter, int num_burn_in
    ) : ctx(ctx), state(init), mem(ctx, seed), factorizer(ctx.config),
        num_iter(num_iter), num_burn_in(num_burn_in) {
            build_spde_matrix(ctx,mem,state.spde_hyperparameters,mem.Q_curr);
            factorizer.submit_matrix(*mem.Q_solver,mem.Q_curr);
            factorizer.factorize_parallel();
        }

    void add_step(std::unique_ptr<Step> s) { steps.push_back(std::move(s)); }
    void add_collector(std::unique_ptr<Recorder> r) { collectors.push_back(std::move(r)); }
    void add_modifier(std::function<void(const SamplerContext &, const SamplerMemory &, const ModelState &)> m) { modifiers.push_back(std::move(m)); }

    ChainOutput run()
    {
        for (int i = 0; i < num_iter; ++i)
        {
            // 1. Preparing all of the matrices for parallel factorization
            for (auto &s : steps)
                s->batch_matrix(ctx, state, mem, factorizer);
            // 2. parallel barrier: all prebatchable factors, each into its own step's storage
            factorizer.factorize_parallel();
            // 3. Gibbs order, serial: RHS assembly + solve + sample only
            for (auto &s : steps)
                s->execute(ctx, state, mem);
            // 4. Recording all of the outputs
            for (auto &c : collectors)
                c->observe(ctx, state, mem);
            // 5. Applying modifiers to the chain 
            for (auto &m : modifiers)
                m(ctx, mem, state);
            mem.curr_iter++;
        }

        ChainOutput out;
        for (auto &c : collectors)
            c->export_to_output(out);
        return out;
    }

private:
    const SamplerContext &ctx;
    ModelState state;
    SamplerMemory mem;
    Factorizer factorizer;
    std::vector<std::unique_ptr<Step>> steps;
    std::vector<std::unique_ptr<Recorder>> collectors;
    std::vector<std::function<void(const SamplerContext &, const SamplerMemory &, const ModelState &)>> modifiers;
    int num_iter, num_burn_in;
};

class ChainRunner
{
    using VecEngines = std::vector<std::unique_ptr<GibbsEngine>>;
    using Configurator = std::function<void(GibbsEngine &)>;

public:
    ChainRunner(
        const SamplerContext &context_,
        const std::vector<ModelState> &init_states,
        uint64_t base_seed,
        size_t num_iter_,
        size_t num_burn_in_)
        : context(context_), num_iter(num_iter_),num_burn_in(num_burn_in_), 
        num_chains(context.config.num_chains),num_threads(context.config.num_chains){
            engines.reserve(num_chains);
            // Constructing the engines
            const uint32_t lo = static_cast<uint32_t>(base_seed);
            const uint32_t hi = static_cast<uint32_t>(base_seed >> 32);

            for(size_t i = 0; i<num_chains;i++){
                std::seed_seq seed({lo,hi,static_cast<uint32_t>(i)});
                engines.push_back(
                    std::make_unique<GibbsEngine>(context,init_states[i],seed,num_iter,num_burn_in)
                );
            }
        }

    void configure(const Configurator &config){
        for (const auto& engine_ptr : engines)
            config(*engine_ptr);
    }
    std::vector<ChainOutput> run(){
        
        std::vector<ChainOutput> outs(num_chains);

        #pragma omp parallel for num_threads(num_threads) if (num_threads > 1) schedule(static)
        for(size_t i = 0;i<num_chains;i++)
            outs[i] = engines[i]->run();
        
        return outs;
    }       

    size_t get_num_iterations() const { return num_iter; }
    size_t get_num_burn_in() const { return num_burn_in; }
    size_t get_num_chains() const { return engines.size(); }

private:
    const SamplerContext &context;
    const size_t num_iter;
    const size_t num_burn_in;
    const size_t num_chains;
    const size_t num_threads;
    VecEngines engines;


};
