#include "FHE_LWE.hpp"

// constructor
FHE_LWE::FHE_LWE() : nttQ1_1(DEGREE1, MODULUS_Q1_1), nttQ1_2(DEGREE1, MODULUS_Q1_2), nttQ2(DEGREE2, MODULUS_Q2), nttQ2Deg1(DEGREE1, MODULUS_Q2)
{
    secretKey_LWE1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));
    secretKey_LWE1_1 =static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));
    secretKey_LWE1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));
    temp_prod_LWE = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE1 * sizeof(uint64_t)));

    secretKey1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    secretKey1_1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    secretKey1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    error1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    error1_1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    error1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    temp_message1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    temp_message1_1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    temp_message1_2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    temp_prod1 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));

    secretKey_LWE2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE_LWE2 * sizeof(uint64_t)));
    secretKey2 = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));
    secretKey2_ntt = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE2 * sizeof(uint64_t)));

    gadget_poly1.resize(GADGET_L_MAX);
    gadget_poly2.resize(GADGET_L_MAX);
    for (int i = 0; i < GADGET_L_MAX; i++)
    {
        gadget_poly1[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
        gadget_poly2[i] = static_cast<uint64_t*>(std::aligned_alloc(64, DEGREE1 * sizeof(uint64_t)));
    }

    KGen_1(secretKey_LWE1, DEGREE_LWE1);
    intel::hexl::EltwiseReduceMod(secretKey_LWE1_1, secretKey_LWE1, DEGREE_LWE1, MODULUS_Q1_1, 1, 1);
    intel::hexl::EltwiseReduceMod(secretKey_LWE1_2, secretKey_LWE1, DEGREE_LWE1, MODULUS_Q1_2, 1, 1);
    KGen_1(secretKey1, DEGREE1);
    intel::hexl::EltwiseReduceMod(secretKey1_1, secretKey1, DEGREE1, MODULUS_Q1_1, 1, 1);
    nttQ1_1.ComputeForward(secretKey1_1, secretKey1_1, 1, 1);
    intel::hexl::EltwiseReduceMod(secretKey1_2, secretKey1, DEGREE1, MODULUS_Q1_1, 1, 1);
    nttQ1_2.ComputeForward(secretKey1_2, secretKey1_2, 1, 1);
    KGen_2(secretKey_LWE2, DEGREE_LWE2);
    KGen_2(secretKey2, DEGREE2);
    nttQ2.ComputeForward(secretKey2_ntt, secretKey2, 1, 1);
}

// destructor
FHE_LWE::~FHE_LWE()
{
    std::free(secretKey_LWE1);
    std::free(secretKey_LWE1_1);
    std::free(secretKey_LWE1_2);
    std::free(temp_prod_LWE);

    std::free(secretKey1);
    std::free(secretKey1_1);
    std::free(secretKey1_2);
    std::free(error1);
    std::free(error1_1);
    std::free(error1_2);
    std::free(temp_message1);
    std::free(temp_message1_1);
    std::free(temp_message1_2);
    std::free(temp_prod1);

    std::free(secretKey_LWE2);
    std::free(secretKey2);
    std::free(secretKey2_ntt);

    for (int i = 0; i < GADGET_L_MAX; i++)
    {
        std::free(gadget_poly1[i]);
        std::free(gadget_poly2[i]);
    }
}

// encodes the message
void FHE_LWE::encode(const int64_t* message, uint64_t* encoded, uint64_t delta ,uint64_t size)
{
    for (size_t i = 0; i < size; ++i)
    {
        encoded[i] = (message[i] < 0) ? static_cast<uint64_t>(MODULUS_Q1 + message[i] * static_cast<int64_t>(delta)) : static_cast<uint64_t>(message[i] * static_cast<int64_t>(delta));
    }
}

// encryption for (n_1, q_1) to LweSample
void FHE_LWE::encryptLweQ1(const uint64_t& message, LweSampleQ1& ciphertext)
{
    // generate random noise
    DiscreteGuassian(0, STDDEVLWE1, MODULUS_Q1, &errorLwe, 1);

    // generate random mask1
    UniformModQ1_1(ciphertext.lweMask1, DEGREE_LWE1);
    intel::hexl::EltwiseMultMod(temp_prod_LWE, ciphertext.lweMask1, secretKey_LWE1_1, DEGREE_LWE1, MODULUS_Q1_1, 1);

    ciphertext.lweVal1 = (sum(temp_prod_LWE, DEGREE_LWE1) + message + errorLwe) % MODULUS_Q1_1;

    // generate random mask2
    UniformModQ1_2(ciphertext.lweMask2, DEGREE_LWE1);
    intel::hexl::EltwiseMultMod(temp_prod_LWE, ciphertext.lweMask2, secretKey_LWE1_2, DEGREE_LWE1, MODULUS_Q1_2, 1);
    ciphertext.lweVal2 = (sum(temp_prod_LWE, DEGREE_LWE1) + message + errorLwe) % MODULUS_Q1_2;
}

// encryption to RLWE over R_{n_1, q_1} 
void FHE_LWE::encryptRlweQ1(const uint64_t* message, RlweSampleQ1& ciphertext)
{
    intel::hexl::EltwiseReduceMod(temp_message1_1, message, DEGREE1, MODULUS_Q1_1, MODULUS_Q1_1, 1);
    nttQ1_1.ComputeForward(temp_message1_1, temp_message1_1, 1, 1);
    intel::hexl::EltwiseReduceMod(temp_message1_2, message, DEGREE1, MODULUS_Q1_2, MODULUS_Q1_2, 1);
    nttQ1_2.ComputeForward(temp_message1_2, temp_message1_2, 1, 1);

    // generate random noise
    DiscreteGuassian(0, STDDEV1, MODULUS_Q1, error1, DEGREE1);
    intel::hexl::EltwiseReduceMod(error1_1, error1, DEGREE1, MODULUS_Q1, MODULUS_Q1_1, 1);
    nttQ1_1.ComputeForward(error1_1, error1_1, 1, 1);
    intel::hexl::EltwiseReduceMod(error1_2, error1, DEGREE1, MODULUS_Q1, MODULUS_Q1_2, 1);
    nttQ1_2.ComputeForward(error1_2, error1_1, 1, 1);

    UniformModQ1_1(ciphertext.rlweMask1, DEGREE1);
    intel::hexl::EltwiseMultMod(temp_prod1, ciphertext.rlweMask1, secretKey1_1, DEGREE1, MODULUS_Q1_1, 1);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal1, temp_prod1, error1_1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal1, ciphertext.rlweVal1, temp_message1_1, DEGREE1, MODULUS_Q1_1); 

    UniformModQ1_2(ciphertext.rlweMask2, DEGREE1);
    intel::hexl::EltwiseMultMod(temp_prod1, ciphertext.rlweMask2, secretKey1_2, DEGREE1, MODULUS_Q1_2, 1);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal2, temp_prod1, error1_2, DEGREE1, MODULUS_Q1_2);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal2, ciphertext.rlweVal2, temp_message1_2, DEGREE1, MODULUS_Q1_2);
}

// encryption to RlevSample over R_{n_1, q_1}
void FHE_LWE::encryptRlevQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t* message, RlevSampleQ1& ciphertext)
{
    lshiftC_Deg1(log_modulus - log_gadget_base * gadget_l, message, temp_message1);
    for (size_t i = 0; i < gadget_l; ++i)
    {
        encryptRlweQ1(temp_message1, ciphertext.rlweSamples[i]);
        lshiftC_I_Deg1(temp_message1, log_gadget_base);
    }
}

// encryption to RgswSample over R_{n_1, q_1}
void FHE_LWE::encryptRgswQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, const uint64_t* message, RgswSampleQ1& ciphertext)
{
    encryptRlevQ1(log_modulus, gadget_l_2, log_gadget_base_2, message, ciphertext.msgRlev);
    nttQ1_1.ComputeForward(temp_message1_1, message, 1, 1);
    intel::hexl::EltwiseMultMod(temp_message1_1, temp_message1_1, secretKey1_1, DEGREE1, MODULUS_Q1_1, 1);
    nttQ1_1.ComputeInverse(temp_message1_1, temp_message1_1, 1, 1);
    nttQ1_2.ComputeForward(temp_message1_2, message, 1, 1);
    intel::hexl::EltwiseMultMod(temp_message1_2, temp_message1_2, secretKey1_2, DEGREE1, MODULUS_Q1_2, 1);
    nttQ1_2.ComputeInverse(temp_message1_2, temp_message1_2, 1, 1);
    CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
    encryptRlevQ1(log_modulus, gadget_l_1, log_gadget_base_1, temp_message1, ciphertext.keyRlev);
}

// half extetnal product for Q1
void FHE_LWE::halfExternalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ1& rlevSample, uint64_t* poly, RlweSampleQ1& result)
{
    GadgetDecompositionDeg1(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
    for (int i = 0; i < gadget_l; i++)
    {
        nttQ1_2.ComputeForward(gadget_poly2[i], gadget_poly1[i], 1, 1);
        nttQ1_1.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
    }
    intel::hexl::EltwiseMultMod(result.rlweMask1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
    intel::hexl::EltwiseMultMod(result.rlweVal1, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);

    intel::hexl::EltwiseMultMod(result.rlweMask2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweMask2, DEGREE1, MODULUS_Q1_1, 2);
    intel::hexl::EltwiseMultMod(result.rlweVal2, gadget_poly2[0], rlevSample.rlweSamples[0].rlweVal2, DEGREE1, MODULUS_Q1_1, 2);

    for (int i = 1; i < gadget_l; i++)
    {
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask1, DEGREE1, MODULUS_Q1_1, 1);
        intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_prod1, DEGREE1, MODULUS_Q1_1);
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal1, DEGREE1, MODULUS_Q1_1, 1);
        intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_prod1, DEGREE1, MODULUS_Q1_1);

        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweMask2, DEGREE1, MODULUS_Q1_1, 2);
        intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_prod1, DEGREE1, MODULUS_Q1_2);
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly2[i], rlevSample.rlweSamples[i].rlweVal2, DEGREE1, MODULUS_Q1_1, 2);
        intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_prod1, DEGREE1, MODULUS_Q1_2);
    }
}

// external product for Q1
void FHE_LWE::externalProdQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result)
{
    nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweMask1, 1, 1);
    nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweMask2, 1, 1);
    CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
    halfExternalProdQ1(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1, result);
    nttQ1_1.ComputeInverse(temp_message1_1, rlweSample.rlweVal1, 1, 1);
    nttQ1_2.ComputeInverse(temp_message1_2, rlweSample.rlweVal2, 1, 1);
    CRT(temp_message1, temp_message1_1, temp_message1_2, DEGREE1);
    halfExternalProdQ1(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1, temp_rlweSampleQ1);
    intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2);
    intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2);
}

// CMux for Q1
void FHE_LWE::CMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample0, RlweSampleQ1& rlweSample1, RlweSampleQ1& result)
{
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ1.rlweMask1, rlweSample1.rlweMask1, rlweSample0.rlweMask1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ1.rlweVal1, rlweSample1.rlweVal1, rlweSample0.rlweVal1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ1.rlweMask2, rlweSample1.rlweMask2, rlweSample0.rlweMask2, DEGREE1, MODULUS_Q1_2);
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ1.rlweVal2, rlweSample1.rlweVal2, rlweSample0.rlweVal2, DEGREE1, MODULUS_Q1_2);
    externalProdQ1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, temp_rlweSampleQ1, result);
    intel::hexl::EltwiseAddMod(result.rlweMask1, result.rlweMask1, rlweSample0.rlweMask1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseAddMod(result.rlweVal1, result.rlweVal1, rlweSample0.rlweVal1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseAddMod(result.rlweMask2, result.rlweMask2, rlweSample0.rlweMask2, DEGREE1, MODULUS_Q1_2);
    intel::hexl::EltwiseAddMod(result.rlweVal2, result.rlweVal2, rlweSample0.rlweVal2, DEGREE1, MODULUS_Q1_2);
}

// DMux for Q1
void FHE_LWE::DMuxQ1(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ1& rgswSample, RlweSampleQ1& rlweSample, RlweSampleQ1& result0, RlweSampleQ1& result1)
{
    externalProdQ1(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, rlweSample, temp_rlweSampleQ1);
    intel::hexl::EltwiseSubMod(result0.rlweMask1, rlweSample.rlweMask1, temp_rlweSampleQ1.rlweMask1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(result0.rlweVal1, rlweSample.rlweVal1, temp_rlweSampleQ1.rlweVal1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(result0.rlweMask2, rlweSample.rlweMask2, temp_rlweSampleQ1.rlweMask2, DEGREE1, MODULUS_Q1_2);
    intel::hexl::EltwiseSubMod(result0.rlweVal2, rlweSample.rlweVal2, temp_rlweSampleQ1.rlweVal2, DEGREE1, MODULUS_Q1_2);

    intel::hexl::EltwiseSubMod(result1.rlweMask1, result0.rlweMask1, rlweSample.rlweMask1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(result1.rlweVal1, result0.rlweVal1, rlweSample.rlweVal1, DEGREE1, MODULUS_Q1_1);
    intel::hexl::EltwiseSubMod(result1.rlweMask2, result0.rlweMask2, rlweSample.rlweMask2, DEGREE1, MODULUS_Q1_2);    
    intel::hexl::EltwiseSubMod(result1.rlweVal2, result0.rlweVal2, rlweSample.rlweVal2, DEGREE1, MODULUS_Q1_2);
}

// encryption to LweSample over (n_2, q_2)
void FHE_LWE::encryptLweQ2(const uint64_t& message, LweSampleQ2& ciphertext)
{
    // generate random noise
    DiscreteGuassian(0, STDDEVLWE2, MODULUS_Q2, &errorLwe, DEGREE_LWE2);

    // generate random mask
    UniformModQ2(ciphertext.lweMask, DEGREE_LWE2);
    intel::hexl::EltwiseMultMod(temp_prod_LWE, ciphertext.lweMask, secretKey_LWE2, DEGREE_LWE2, MODULUS_Q2, 1);

    ciphertext.lweVal = (sum(temp_prod_LWE, DEGREE_LWE2) + message + errorLwe) % MODULUS_Q2;
}

// encryption to RLWE over R_{n_2, q_2}
void FHE_LWE::encryptRlweQ2(const uint64_t* message, RlweSampleQ2& ciphertext)
{
    nttQ2.ComputeForward(temp_message1_1, message, 1, 1);

    // generate random noise
    DiscreteGuassian(0, STDDEV2, MODULUS_Q2, error1, DEGREE2);
    nttQ2.ComputeForward(error1_1, error1_1, 1, 1);

    UniformModQ1_1(ciphertext.rlweMask, DEGREE2);
    intel::hexl::EltwiseMultMod(temp_prod1, ciphertext.rlweMask, secretKey2, DEGREE2, MODULUS_Q2, 1);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal, temp_prod1, error1_1, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseAddMod(ciphertext.rlweVal, ciphertext.rlweVal, temp_message1_1, DEGREE2, MODULUS_Q2); 
}

// encryption to RlevSample over R_{n_2, q_2}
void FHE_LWE::encryptRlevQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, const uint64_t* message, RlevSampleQ2& ciphertext)
{
    lshiftC_Deg2(log_modulus - log_gadget_base * gadget_l, message, temp_message1);
    for (size_t i = 0; i < gadget_l; ++i)
    {
        encryptRlweQ2(temp_message1, ciphertext.rlweSamples[i]);
        lshiftC_I_Deg2(temp_message1, log_gadget_base);
    }
}

// encryption to RgswSample over R_{n_2, q_2}
void FHE_LWE::encryptRgswQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, const uint64_t* message, RgswSampleQ2& ciphertext)
{
    encryptRlevQ2(log_modulus, gadget_l_2, log_gadget_base_2, message, ciphertext.msgRlev);
    nttQ2.ComputeForward(temp_message1_1, message, 1, 1);
    intel::hexl::EltwiseMultMod(temp_message1, temp_message1_1, secretKey2_ntt, DEGREE2, MODULUS_Q2, 1);
    nttQ2.ComputeInverse(temp_message1, temp_message1, 1, 1);
    encryptRlevQ2(log_modulus, gadget_l_1, log_gadget_base_1, temp_message1, ciphertext.keyRlev);
}

// half product for Q2
void FHE_LWE::halfExternalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RlevSampleQ2& rlevSample, uint64_t* poly, RlweSampleQ2& result)
{
    GadgetDecompositionDeg2(poly, gadget_poly1, log_modulus, gadget_l, log_gadget_base);
    for (int i = 0; i < gadget_l; i++)
    {
        nttQ2.ComputeForward(gadget_poly1[i], gadget_poly1[i], 1, 1);
    }
    intel::hexl::EltwiseMultMod(result.rlweMask, gadget_poly1[0], rlevSample.rlweSamples[0].rlweMask, DEGREE2, MODULUS_Q2, 1);
    intel::hexl::EltwiseMultMod(result.rlweVal, gadget_poly1[0], rlevSample.rlweSamples[0].rlweVal, DEGREE2, MODULUS_Q2, 1);

    for (int i = 1; i < gadget_l; i++)
    {
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweMask, DEGREE2, MODULUS_Q2, 1);
        intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_prod1, DEGREE2, MODULUS_Q2);
        intel::hexl::EltwiseMultMod(temp_prod1, gadget_poly1[i], rlevSample.rlweSamples[i].rlweVal, DEGREE2, MODULUS_Q2, 1);
        intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_prod1, DEGREE2, MODULUS_Q2);
    }
}

// external product for Q2
void FHE_LWE::externalProdQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result)
{
    nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweMask, 1, 1);
    halfExternalProdQ2(log_modulus, gadget_l_1, log_gadget_base_1, rgswSample.keyRlev, temp_message1_1, result);
    nttQ2.ComputeInverse(temp_message1_1, rlweSample.rlweVal, 1, 1);
    halfExternalProdQ2(log_modulus, gadget_l_2, log_gadget_base_2, rgswSample.msgRlev, temp_message1_1, temp_rlweSampleQ2);
    intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, temp_rlweSampleQ2.rlweMask, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, temp_rlweSampleQ2.rlweVal, DEGREE2, MODULUS_Q2);
}

// CMux for Q2
void FHE_LWE::CMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample0, RlweSampleQ2& rlweSample1, RlweSampleQ2& result)
{
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweMask, rlweSample1.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseSubMod(temp_rlweSampleQ2.rlweVal, rlweSample1.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
    externalProdQ2(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, temp_rlweSampleQ2, result);
    intel::hexl::EltwiseAddMod(result.rlweMask, result.rlweMask, rlweSample0.rlweMask, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseAddMod(result.rlweVal, result.rlweVal, rlweSample0.rlweVal, DEGREE2, MODULUS_Q2);
}

 // DMux for Q2
void FHE_LWE::DMuxQ2(const uint64_t log_modulus, const uint64_t gadget_l_1, const uint64_t gadget_l_2, const uint64_t log_gadget_base_1, const uint64_t log_gadget_base_2, RgswSampleQ2& rgswSample, RlweSampleQ2& rlweSample, RlweSampleQ2& result0, RlweSampleQ2& result1)
{
    externalProdQ2(log_modulus, gadget_l_1, gadget_l_2, log_gadget_base_1, log_gadget_base_2, rgswSample, rlweSample, temp_rlweSampleQ2);
    intel::hexl::EltwiseSubMod(result0.rlweMask, rlweSample.rlweMask, temp_rlweSampleQ2.rlweMask, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseSubMod(result0.rlweVal, rlweSample.rlweVal, temp_rlweSampleQ2.rlweVal, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseSubMod(result1.rlweMask, result0.rlweMask, rlweSample.rlweMask, DEGREE2, MODULUS_Q2);
    intel::hexl::EltwiseSubMod(result1.rlweVal, result0.rlweVal, rlweSample.rlweVal, DEGREE2, MODULUS_Q2);
}

// generate ring switching key
void FHE_LWE::GenRsk(const uint64_t log_modulus, const uint64_t gadget_l, const uint64_t log_gadget_base, RingSwitchingKey& rsk)
{
    ring_embed(DEGREE2, DEGREE1, 0, secretKey2, temp_message1);
    nttQ2.ComputeForward(temp_message1, temp_message1, 1, 1);
    lshiftC_Deg1(log_modulus - log_gadget_base * gadget_l, temp_message1_2, secretKey1);
    for (size_t i = 0; i < gadget_l; ++i)
    {
        nttQ2Deg1.ComputeForward(temp_message1_1, temp_message1_2, 1, 1);

        // generate random noise
        DiscreteGuassian(0, STDDEV2, MODULUS_Q2, error1, DEGREE1);
        nttQ2Deg1.ComputeForward(error1_1, error1_1, 1, 1);

        UniformModQ1_1(rsk.rlweMasks[i], DEGREE1);
        intel::hexl::EltwiseMultMod(temp_prod1, rsk.rlweMasks[i], temp_message1, DEGREE1, MODULUS_Q2, 1);
        intel::hexl::EltwiseAddMod(rsk.rlweVals[i], temp_prod1, error1_1, DEGREE1, MODULUS_Q2);
        intel::hexl::EltwiseAddMod(rsk.rlweVals[i], rsk.rlweVals[i], temp_message1_1, DEGREE1, MODULUS_Q2); 
        lshiftC_I_Deg2(temp_message1_2, log_gadget_base);
    }
}