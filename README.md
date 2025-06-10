# VIA: Communication-Efficient Single Server Private Information Retrieval

**This repository implements VIA and VIA-C**, two communication-efficient single-server Private Information Retrieval (PIR) schemes proposed in the paper _"VIA: Communication-Efficient Single Server Private Information Retrieval"_.

**Warning: This code is intended for academic research purposes only.**

## Running Experiments

Our code was tested on an Ubuntu 22.04.5 workstation with an AMD Ryzen 9 9950X CPU and 128GB RAM. The code was compiled using clang 19.1.7 with support for AVX-512DQ and AVX-512IFMA52 instruction sets.

### Prerequisites
1. **Build Intel HEXL Library**: Install and build the [Intel HEXL library](https://github.com/IntelLabs/hexl) before compiling this code. Place the resulting static library `libhexl.a` in the repository's `lib` directory.
2. **Verify CPU Support**: Ensure your CPU supports the AVX-512DQ instruction set.
3. **Enable AVX-512**: Enable AVX-512 support in your compiler settings during compilation.

### Configuration
Adjust parameters by modifying `core.hpp` to test PIR under different configurations.

### Performance Tests

The performance evaluation offers the following test variants:  
- `VIA`  
- `VIA with Blinded Extraction`  
- `VIA-C`  
- `VIA-C with Blinded Extraction`  

**Execution Notes:**  When running tests:  
1. Select either VIA or VIA-C  
2. Choose whether to enable Blinded Extraction  
3. Specify the number of test iterations  

#### Memory Handling for Large Databases
Due to memory constraints:  
- **VIA**: For databases > 4GB, only 4GB database of memory will be allocated. The First Dimension operations will simulate memory usage for the full database.  
- **VIA-C**: For databases >2GB, only 2GB database of memory will be allocated. First Dimension operations will simulate memory usage for the full database.  



