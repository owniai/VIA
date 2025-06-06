#pragma once
#include "core.hpp"

class LweSampleQ1
{
public:
    // constructor
    LweSampleQ1()
    {
        lweMask1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));
        lweMask2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));
    }

    // destructor
    ~LweSampleQ1()
    {
        // std::free(lweMask1);
        // std::free(lweMask2);
    }

    // variables
    uint64_t* lweMask1;
    uint64_t lweVal1;
    uint64_t* lweMask2;
    uint64_t lweVal2;
};

class LweSampleQ2
{
public:
    // constructor
    LweSampleQ2()
    {
        lweMask = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE2 * sizeof(uint64_t)));
    }

    // destructor
    ~LweSampleQ2()
    {
        // std::free(lweMask);
    }

    // variables
    uint64_t* lweMask;
    uint64_t lweVal;
};