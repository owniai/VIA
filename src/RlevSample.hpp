#pragma once
#include "RlweSample.hpp"

class RlevSampleQ1
{
public:
    // Constructor
    RlevSampleQ1(uint64_t l)
    {
        rlweSamples.resize(l);
    }

    // Destructor
    ~RlevSampleQ1() = default;

    // variables
    std::vector<RlweSampleQ1> rlweSamples;
};

class RlevSampleQ2Deg1
{
public:
    // Constructor
    RlevSampleQ2Deg1(uint64_t l)
    {
        rlweSamples.resize(l);
    }

    // Destructor
    ~RlevSampleQ2Deg1() = default;

    // variables
    std::vector<RlweSampleQ2Deg1> rlweSamples;
};

class RlevSampleQ2
{
public:
    // Constructor
    RlevSampleQ2(uint64_t l)
    {
        rlweSamples.resize(l);
    }

    // Destructor
    ~RlevSampleQ2() = default;

    // variables
    std::vector<RlweSampleQ2> rlweSamples;
};

class RingSwitchingKey
{
public:
    // Constructor
    RingSwitchingKey(uint64_t l)
    {
        rlweMasks.resize(l);
        rlweVals.resize(l);
        for (uint64_t i = 0; i < l; i++)
        {
            rlweMasks[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
            rlweVals[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        }
    }

    // Destructor
    ~RingSwitchingKey()
    {
        for (uint64_t i = 0; i < rlweMasks.size(); i++)
        {
            std::free(rlweMasks[i]);
            std::free(rlweVals[i]);
        }
    }

    // variables
    std::vector<uint64_t*> rlweMasks;
    std::vector<uint64_t*> rlweVals;
};

