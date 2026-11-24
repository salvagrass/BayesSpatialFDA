#pragma once

#include "SamplingUtils.hpp"
#include "GibbsSteps.hpp"
#include "GibbsHooks.hpp"

class GibbsEngine{
public:
    ChainOutput run();


private:
    std::vector<Step> steps;
    std::vector<Recorder> recorders;
};



class ChainRunner{
    using VecEngines = std::vector<std::unique_ptr<GibbsEngine>>;
    using HookConfigurator = std::function<void(GibbsEngine&)>;

public:
    ChainRunner(
        const SamplerContext& context_,
        const std::vector<ModelState>& init_states, 
        uint64_t base_seed,
        size_t num_iter_,
        size_t num_burn_in_
    );


    void configure_hook(const HookConfigurator& config);
    std::vector<ChainOutput> run();


    size_t get_num_iterations() const { return num_iter; }
    size_t get_num_burn_in() const    { return num_burn_in; }
    size_t get_num_chains() const     { return engines.size(); }

private:
    const SamplerContext& context;
    const size_t num_iter;
    const size_t num_burn_in;
    const size_t num_chains;
    VecEngines engines;
};



