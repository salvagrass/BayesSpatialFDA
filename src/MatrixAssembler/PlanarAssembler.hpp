#pragma once

#include "SPDEAssembler.hpp"


class PlanarAssembler : public SPDEAssembler{
    public:
        using SPDEAssembler::SPDEAssembler;

        void assemble_matrices() override;
};