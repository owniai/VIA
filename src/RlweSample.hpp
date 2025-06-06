#pragma once
#include "core.hpp"

class RlweSampleQ1
{
public:
	// constructor
	RlweSampleQ1()
	{
		rlweMask1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
		rlweVal1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
		rlweMask2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
		rlweVal2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
	}
	// destructor
	~RlweSampleQ1() = default;

	// variables
    uint64_t* rlweMask1;
	uint64_t* rlweVal1;
	uint64_t* rlweMask2;
	uint64_t* rlweVal2;
};

class RlweSampleQ2Deg1
{
public:
	// constructor
	RlweSampleQ2Deg1()
	{
		rlweMask = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
		rlweVal = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
	}

	// destructor
	~RlweSampleQ2Deg1() = default;

	// variables
    uint64_t* rlweMask;
	uint64_t* rlweVal;
};

class RlweSampleQ2
{
public:
	// constructor
	RlweSampleQ2()
	{
		rlweMask = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
		rlweVal = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
	}

	// destructor
	~RlweSampleQ2() = default;

	// variables
    uint64_t* rlweMask;
	uint64_t* rlweVal;
};