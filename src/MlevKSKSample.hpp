#pragma once
#include "core.hpp"
#include "MlweSample.hpp"

class MlevKSKSampleQ1
{
public:
    MlevKSKSampleQ1(uint64_t k, uint64_t deg, uint64_t l)
    {
        for (uint64_t i = 0; i < l*k; i++)
        {
            mlweSamples.emplace_back(MlweSampleQ1(k, deg));
        }
    }

    ~MlevKSKSampleQ1() = default;

    std::vector<MlweSampleQ1> mlweSamples;
};

class MlevKSKSampleQ2
{
public:
    MlevKSKSampleQ2(uint64_t k, uint64_t deg, uint64_t l)
    {
        for (uint64_t i = 0; i < l*k; i++)
        {
            mlweSamples.emplace_back(MlweSampleQ2(k, deg));
        }
    }

    ~MlevKSKSampleQ2() = default;

    std::vector<MlweSampleQ2> mlweSamples;
};