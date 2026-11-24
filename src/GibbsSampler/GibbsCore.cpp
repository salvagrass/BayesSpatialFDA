#include "GibbsCore.hpp"
#include <stdexcept>


ChainRunner::ChainRunner(
    const SamplerContext& context_, 
    const std::vector<ModelState>& init_states,
    uint64_t base_seed,
    size_t num_iter_,
    size_t num_burn_in_
): context(context_), num_iter(num_iter_),num_burn_in(num_burn_in_), num_chains(context.config.num_chains){

    // Error Checking
    if (init_states.size() != num_chains) {
        throw std::invalid_argument("init_states size must match num_chains.");
    }
    if (num_burn_in >= num_iter) {
        throw std::invalid_argument("num_burn_in must be strictly less than num_iter.");
    }

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


void ChainRunner::configure_hook(const HookConfigurator& config){
    for(const auto& engine_ptr : engines)
        config(*engine_ptr);
}


std::vector<ChainOutput> ChainRunner::run(){
    std::vector<ChainOutput> outs(num_chains);

    for(size_t idx=0; idx<num_chains;idx++)
        outs[idx] = engines[idx]->run(); 

    return outs;
}