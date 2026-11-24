#pragma once

#include "SamplingUtils.hpp"
#include "GibbsSteps.hpp"
#include "GibbsHooks.hpp"

class GibbsEngine
{
public:
    GibbsEngine(const SamplerContext &ctx, const ModelState &init,
                const SamplerMemory &mem,
                std::seed_seq &seed, int num_iter, int num_burn_in,
                int factor_threads = 1)
        : ctx_(ctx), state_(init), mem_(ctx, seed),
          factorizer_(factor_threads),
          num_iter_(num_iter), num_burn_in_(num_burn_in) {}

    void add_step(std::unique_ptr<Step> s) { steps_.push_back(std::move(s)); }
    void add_collector(std::unique_ptr<Recorder> r) { collectors_.push_back(std::move(r)); }
    void add_collector(std::function<void(const SamplerContext &, const SamplerMemory &, const ModelState &)> m) { modifiers.push_back(std::move(m)); }

    ChainOutput run()
    {

        for (int i = 0; i < num_iter_; ++i)
        {
            for (auto &s : steps_)
                s->batch_matrix(ctx_, state_, mem_, factorizer_);

            // 2. parallel barrier: all prebatchable factors, each into its own step's storage
            factorizer_.factorize_parallel();

            // 3. Gibbs order, serial: RHS assembly + solve + sample only
            for (auto &s : steps_)
                s->execute(ctx_, state_, mem_);

            // 5. thinned, post-burn-in recording
            for (auto &c : collectors_)
                c->observe(ctx_, state_, mem_, i);

            //
            for (auto &m : modifiers)
                m(ctx_, mem_, state_);
        }

        ChainOutput out;
        for (auto &c : collectors_)
            c->export_to_output(out);
        return out;
    }

private:
    const SamplerContext &ctx_;
    ModelState state_;
    SamplerMemory mem_;
    Factorizer factorizer_;
    std::vector<std::unique_ptr<Step>> steps_;
    std::vector<std::unique_ptr<Recorder>> collectors_;
    std::vector<std::function<void(const SamplerContext &, const SamplerMemory &, const ModelState &)>> modifiers;
    int num_iter_, num_burn_in_;
};

class ChainRunner
{
    using VecEngines = std::vector<std::unique_ptr<GibbsEngine>>;
    using HookConfigurator = std::function<void(GibbsEngine &)>;

public:
    ChainRunner(
        const SamplerContext &context_,
        const std::vector<ModelState> &init_states,
        uint64_t base_seed,
        size_t num_iter_,
        size_t num_burn_in_);

    void configure_hook(const HookConfigurator &config);
    std::vector<ChainOutput> run();

    size_t get_num_iterations() const { return num_iter; }
    size_t get_num_burn_in() const { return num_burn_in; }
    size_t get_num_chains() const { return engines.size(); }

private:
    const SamplerContext &context;
    const size_t num_iter;
    const size_t num_burn_in;
    const size_t num_chains;
    VecEngines engines;
};
