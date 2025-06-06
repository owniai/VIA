#pragma once
#include "core.hpp"
#include "LweSample.hpp"
#include "RgswSample.hpp"
#include "functions.hpp"

class FHE_LWE
{
public:

    // constructor
    FHE_LWE();

    // destructor
    ~FHE_LWE();

    // encodes the message
    void encode(const int64_t* message, uint64_t* encoded, uint64_t delta, uint64_t size);

    // encryption to LweSample over (n_1, q_1)
    void encryptLweQ1(const uint64_t& message, LweSampleQ1& ciphertext);

    // encryption to RLWE over R_{n_1, q_1} 
    void encryptRlweQ1(const uint64_t* message, RlweSampleQ1& ciphertext);

    // encryption to RlevSample over R_{n_1, q_1}
    void encryptRlevQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t* message, RlevSampleQ1& ciphertext);

    // encryption to RgswSample over R_{n_1, q_1}
    void encryptRgswQ1(const uint64_t log_modulus, const uint64_t gadget_l_A, const uint64_t gadget_l_B, const uint64_t log_gadget_base_A, const uint64_t log_gadget_base_B, const uint64_t* message, RgswSampleQ1& ciphertext);

    // half external product for Q1
    void halfExternalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ1& rlevSample, uint64_t* poly, RlweSampleQ1& result);

    // external product for Q1
    void externalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result);

    // CMux for Q1
    void CMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample0, RlweSampleQ1& rlweSample1, RlweSampleQ1& result);

    // DMux for Q1
    void DMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result0, RlweSampleQ1& result1);

    // encryption to LweSample over (n_2, q_3)
    void encryptLweQ2(const uint64_t& message, LweSampleQ2& ciphertext);

    // encryption to RLWE over R_{n_2, q_3}
    void encryptRlweQ2(const uint64_t* message, RlweSampleQ2& ciphertext);

    // encryption to RlevSample over R_{n_2, q_3}
    void encryptRlevQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t* message, RlevSampleQ2& ciphertext);

    // encryption to RgswSample over R_{n_2, q_3}
    void encryptRgswQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, const uint64_t* message, RgswSampleQ2& ciphertext);

    // half product for Q2
    void halfExternalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ2& rlevSample, uint64_t* poly, RlweSampleQ2& result);

    // external product for Q2
    void externalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result);

    // CMux for Q2
    void CMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample0, RlweSampleQ2& rlweSample1, RlweSampleQ2& result);

    // DMux for Q2
    void DMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result0, RlweSampleQ2& result1);

    // generate ring switching key
    void GenRsk(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RingSwitchingKey& rsk);

    // variables
    intel::hexl::NTT nttQ1_1;
    intel::hexl::NTT nttQ1_2;
    intel::hexl::NTT nttQ2;
    intel::hexl::NTT nttQ2Deg1;

    uint64_t* secretKey_LWE1;
    uint64_t* secretKey_LWE1_1;
    uint64_t* secretKey_LWE1_2;
    uint64_t errorLwe;
    uint64_t decRes1;
    uint64_t decRes2;
    int64_t decRes;
    uint64_t* temp_prod_LWE;

    uint64_t* secretKey1;
    uint64_t* secretKey1_1;
    uint64_t* secretKey1_2;

    uint64_t* error1;
    uint64_t* error1_1;
    uint64_t* error1_2;
    uint64_t* temp_message1;
    uint64_t* temp_message1_1;
    uint64_t* temp_message1_2;
    uint64_t* temp_prod1;

    uint64_t* secretKey_LWE2;
    uint64_t* secretKey2;
    uint64_t* secretKey2_ntt;

    std::vector<uint64_t*> gadget_poly1;
    std::vector<uint64_t*> gadget_poly2;
    RlweSampleQ1 temp_rlweSampleQ1;
    RlweSampleQ2 temp_rlweSampleQ2;
};