#pragma once
#include "RlevSample.hpp"

class RgswSampleQ1
{
public:
    // constructor
    RgswSampleQ1(uint64_t l_1, uint64_t l_2) : msgRlev(l_2), keyRlev(l_1) {}

    // destructor
    ~RgswSampleQ1() = default;

    // variables
    RlevSampleQ1 msgRlev;
    RlevSampleQ1 keyRlev;
};

class RgswSampleQ2Deg1
{
public:
    // constructor
    RgswSampleQ2Deg1(uint64_t l_1, uint64_t l_2) : msgRlev(l_2), keyRlev(l_1) {}

    // destructor
    ~RgswSampleQ2Deg1() = default;

    // variables
    RlevSampleQ2Deg1 msgRlev;
    RlevSampleQ2Deg1 keyRlev;
};

class RgswSampleQ2
{
public:
    // constructor
    RgswSampleQ2(uint64_t l_1, uint64_t l_2) : msgRlev(l_2), keyRlev(l_1) {}

    // destructor
    ~RgswSampleQ2() = default;
    
    // variables
    RlevSampleQ2 msgRlev;
    RlevSampleQ2 keyRlev;
};