#pragma once
#include "core.hpp"

class MlweSampleQ1
{
public:
    MlweSampleQ1(uint64_t k, uint64_t deg)
    {
        mlweMasks1.resize(k);
        mlweMasks2.resize(k);
        for (uint64_t i = 0; i < k; i++)
        {
            mlweMasks1[i] = static_cast<uint64_t*>(std::aligned_alloc(64, deg * sizeof(uint64_t)));
            mlweMasks2[i] = static_cast<uint64_t*>(std::aligned_alloc(64, deg * sizeof(uint64_t)));
        }
        mlweVal1 = static_cast<uint64_t*>(std::aligned_alloc(64, sizeof(uint64_t)));
        mlweVal2 = static_cast<uint64_t*>(std::aligned_alloc(64, sizeof(uint64_t)));
    }

    ~MlweSampleQ1()
    {
        for (uint64_t i = 0; i < mlweMasks1.size(); i++)
        {
            std::free(mlweMasks1[i]);
            std::free(mlweMasks2[i]);
        }
        std::free(mlweVal1);
        std::free(mlweVal2);
    }

    std::vector<uint64_t*> mlweMasks1;
    std::vector<uint64_t*> mlweMasks2;
    uint64_t* mlweVal1;
    uint64_t* mlweVal2;
};

class MlweSampleQ2
{
public:
    MlweSampleQ2(uint64_t k, uint64_t deg)
    {
        mlweMasks.resize(k);
        for (uint64_t i = 0; i < k; i++)
        {
            mlweMasks[i] = static_cast<uint64_t*>(std::aligned_alloc(64, deg * sizeof(uint64_t)));
        }
        mlweVal = static_cast<uint64_t*>(std::aligned_alloc(64, sizeof(uint64_t)));
    }

    ~MlweSampleQ2()
    {
        for (uint64_t i = 0; i < mlweMasks.size(); i++)
        {
            std::free(mlweMasks[i]);
        }
        std::free(mlweVal);
    }

    std::vector<uint64_t*> mlweMasks;
    uint64_t* mlweVal;
};