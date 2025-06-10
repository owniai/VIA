#pragma once
#include "core.hpp"
#include <random>
#include <immintrin.h>

// set seeds sequence
constexpr static const std::initializer_list<uint32_t> seeds = { 142, 258, 1598, 0711, 0212, 2016, 1536 };
static std::seed_seq seq(seeds);

// random number generator
static std::mt19937_64 gen(seq);

// normal distribution with mean 0 and stddev 1
static std::normal_distribution<float> normal_dist(0.0f, 1.0f);
// secret key generator 1 distribution
static std::uniform_int_distribution<int64_t> KGen_1_dist(-HALF_MODULUS_KGEN1, HALF_MODULUS_KGEN1);
// secret key generator 2 distribution
static std::uniform_int_distribution<int64_t> KGen_2_dist(-HALF_MODULUS_KGEN2, HALF_MODULUS_KGEN2);
// secret key generator 1C distribution
static std::uniform_int_distribution<int64_t> KGen_1C_dist(-HALF_MODULUS_KGEN1C, HALF_MODULUS_KGEN1C);
// secret key generator 2C distribution
static std::uniform_int_distribution<int64_t> KGen_2C_dist(-HALF_MODULUS_KGEN2C, HALF_MODULUS_KGEN2C);
// uniform distribution mod Q1_1
static std::uniform_int_distribution<int64_t> UniformModQ1_1_dist(-HALF_MODULUS_Q1_1, HALF_MODULUS_Q1_1);
// uniform distribution mod Q1_2
static std::uniform_int_distribution<int64_t> UniformModQ1_2_dist(-HALF_MODULUS_Q1_2, HALF_MODULUS_Q1_2);
// Uniform distribution mod Q2
static std::uniform_int_distribution<int64_t> UniformModQ2_dist(-HALF_MODULUS_Q2, HALF_MODULUS_Q2);
// Uniform distribution mod P
static std::uniform_int_distribution<int64_t> UniformModP_dist(-HALF_MODULUS_P, HALF_MODULUS_P-1);
// uniform distribution mod Q1_1C
static std::uniform_int_distribution<int64_t> UniformModQ1_1C_dist(-HALF_MODULUS_Q1_1C, HALF_MODULUS_Q1_1C);
// uniform distribution mod Q1_2C
static std::uniform_int_distribution<int64_t> UniformModQ1_2C_dist(-HALF_MODULUS_Q1_2C, HALF_MODULUS_Q1_2C);
// Uniform distribution mod Q2C
static std::uniform_int_distribution<int64_t> UniformModQ2C_dist(-HALF_MODULUS_Q2C, HALF_MODULUS_Q2C);

// Discrete Gaussian Distribution mod modulus
void DiscreteGuassian(float mean, float stddev, uint64_t modulus, uint64_t* result, uint64_t size);

// secret key generator 1
void KGen_1(uint64_t* result, uint64_t size);

// secret key generator 2
void KGen_2(uint64_t* result, uint64_t size);

// secret key generator 1
void KGen_1C(uint64_t* result, uint64_t size);

// secret key generator 2
void KGen_2C(uint64_t* result, uint64_t size);

//  uniform mod Q1_1
void UniformModQ1_1(uint64_t* result, uint64_t size);

//  uniform mod Q1_2
void UniformModQ1_2(uint64_t* result, uint64_t size);

//  uniform mod Q2
void UniformModQ2(uint64_t* result, uint64_t size);

//  uniform mod P
void UniformModP(uint64_t* result, uint64_t size);

//  uniform mod Q1_1C
void UniformModQ1_1C(uint64_t* result, uint64_t size);

//  uniform mod Q1_2C
void UniformModQ1_2C(uint64_t* result, uint64_t size);

//  uniform mod Q2C
void UniformModQ2C(uint64_t* result, uint64_t size);

// the sum of elements in Src
uint64_t sum(const uint64_t* Src, uint64_t size);

// solve CRT equations
void CRT(uint64_t& result, const uint64_t& op1, const uint64_t& op2);
void CRT(uint64_t* result, uint64_t* op1, uint64_t* op2, uint64_t size);

// solve CRT equations for VIA-C
void CRT_C(uint64_t& result, const uint64_t& op1, const uint64_t& op2);
void CRT_C(uint64_t* result, uint64_t* op1, uint64_t* op2, uint64_t size);

// Gadget Decomposition for poly. of degree1
void GadgetDecompositionDeg1(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus ,const uint64_t gadget_l, const uint64_t log_gadget_base);

// Gadget Decomposition for poly. of degree2
void GadgetDecompositionDeg2(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus ,const uint64_t gadget_l, const uint64_t log_gadget_base);

// Gadget Decomposition for poly.
void GadgetDecompositionDegAll(uint64_t* poly, std::vector<uint64_t*> result, const uint64_t log_modulus ,const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t deg);

// ring embedding
void ring_embed(const uint64_t deg_src, const uint64_t deg_dst, const uint64_t idx, const uint64_t* src, uint64_t* dst);

// ring projection
void ring_project(const uint64_t idx, const uint64_t* src, uint64_t* dst);

// modulus switching
void modSwitch(const uint64_t modulus_src, const uint64_t modulus_dst, uint64_t* Src, uint64_t* Dst, const uint64_t size);

// printInfoVIA
void printInfoVIA(uint64_t IsBlindedExtraction);

// printInfoVIA
void printInfoVIA_C(uint64_t IsBlindedExtraction);

// AVX512 multiply Constant
void mulC_Deg1(const uint64_t val, uint64_t* Src, uint64_t* Dst);

// AVX512 In-place multiply Constant
void mulC_I_Deg1(const uint64_t val, uint64_t* SrcDst);

// AVX512 addition
void add_Deg1(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst);

// AVX512 In-place addition
void add_I_Deg1(const uint64_t* Src, uint64_t* SrcDst);

// AVX512 subtraction
void sub_Deg1(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst);

// AVX512 In-place subtraction
void sub_I_Deg1(uint64_t* SrcDst, const uint64_t* Src);

// AVX512 negative
void negative_Deg1(uint64_t* Src, uint64_t* Dst);

// AVX512 In-place negative
void negative_I_Deg1(uint64_t* SrcDst);

// AVX512 copy
void copy_Deg1(const uint64_t* Src, uint64_t* Dst);

// AVX512 add Constant
void addC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place add Constant
void addC_I_Deg1(const uint64_t val, uint64_t* Dst);

// AVX512 shift right
void rshiftC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place shift right
void rshiftC_I_Deg1(uint64_t* SrcDst, uint64_t shift);

// AVX512 In-place shift right signed
void rshiftC_Signed_I_Deg1(int64_t* SrcDst, uint64_t shift);

// AVX512 shift left
void lshiftC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place shift left
void lshiftC_I_Deg1(uint64_t* SrcDst, uint64_t shift);

// AVX512 logic and constant
void andC_Deg1(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 multiply Constant
void mulC_Deg2(const uint64_t val, uint64_t* Src, uint64_t* Dst);

// AVX512 In-place multiply Constant
void mulC_I_Deg2(const uint64_t val, uint64_t* SrcDst);

// AVX512 addition
void add_Deg2(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst);

// AVX512 In-place addition
void add_I_Deg2(const uint64_t* Src, uint64_t* SrcDst);

// AVX512 subtraction
void sub_Deg2(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst);

// AVX512 In-place subtraction
void sub_I_Deg2(uint64_t* SrcDst, const uint64_t* Src);

// AVX512 negative
void negative_Deg2(uint64_t* Src, uint64_t* Dst);

// AVX512 In-place negative
void negative_I_Deg2(uint64_t* SrcDst);

// AVX512 copy
void copy_Deg2(const uint64_t* Src, uint64_t* Dst);

// AVX512 add Constant
void addC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place add Constant
void addC_I_Deg2(const uint64_t val, uint64_t* Dst);

// AVX512 shift right
void rshiftC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place shift right
void rshiftC_I_Deg2(uint64_t* SrcDst, uint64_t shift);

// AVX512 In-place shift right signed
void rshiftC_Signed_I_Deg2(int64_t* SrcDst, uint64_t shift);

// AVX512 shift left
void lshiftC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// AVX512 In-place shift left
void lshiftC_I_Deg2(uint64_t* SrcDst, uint64_t shift);

// AVX512 logic and constant
void andC_Deg2(const uint64_t val, const uint64_t* Src, uint64_t* Dst);

// -----------------
// AVX512 multiply Constant
void mulC_DegAll(const uint64_t val, uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 In-place multiply Constant
void mulC_I_DegAll(const uint64_t val, uint64_t* SrcDst, const uint64_t size);

// AVX512 addition
void add_DegAll(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst, const uint64_t size);

// AVX512 In-place addition
void add_I_DegAll(const uint64_t* Src, uint64_t* SrcDst, const uint64_t size);

// AVX512 subtraction
void sub_DegAll(const uint64_t* Src1, const uint64_t* Src2, uint64_t* Dst, const uint64_t size);

// AVX512 In-place subtraction
void sub_I_DegAll(uint64_t* SrcDst, const uint64_t* Src, const uint64_t size);

// AVX512 negative
void negative_DegAll(uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 In-place negative
void negative_I_DegAll(uint64_t* SrcDst, const uint64_t size);

// AVX512 copy
void copy_DegAll(const uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 add Constant
void addC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 In-place add Constant
void addC_I_DegAll(const uint64_t val, uint64_t* Dst, const uint64_t size);

// AVX512 shift right
void rshiftC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 In-place shift right
void rshiftC_I_DegAll(uint64_t* SrcDst, uint64_t shift, const uint64_t size);

// AVX512 In-place shift right signed
void rshiftC_Signed_I_DegAll(int64_t* SrcDst, uint64_t shift, const uint64_t size);

// AVX512 shift left
void lshiftC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size);

// AVX512 In-place shift left
void lshiftC_I_DegAll(uint64_t* SrcDst, uint64_t shift, const uint64_t size);

// AVX512 logic and constant
void andC_DegAll(const uint64_t val, const uint64_t* Src, uint64_t* Dst, const uint64_t size);
