
# **BFNET: A Modern Fortran Machine Learning Potential Code**

BFNET (**Behler-Parrinello neural NETwork potential**) is an open-source, modern Fortran implementation of machine learning potential construction. It is inspired by [aenet](http://ann.atomistic.net) and aims to provide a more maintainable and extensible alternative while introducing new features and optimizations.


**Note:** BFNET is still under development and does not yet support all features of aenet. However, it is designed to **extend beyond aenet** by incorporating new capabilities.

**Architecture**

BFNET is structured into **three key components**, leveraging **modern Fortran** for high-performance computation and **Julia** (under development) for flexible model training:

1. **Descriptor Generation (Modern Fortran)**

- Computes atomic environment descriptors based on Behler-Parrinello-type symmetry functions. Now only Chebyshev descriptor (compatible with aenet v2.0.4 not 2.0.3) is supported.

2. **Neural Network Training (Julia, in development)**

- The training module is currently under development and will be released in a separate repository.

- In the meantime, **aenet can be used for training**, as BFNET’s descriptor format is compatible with aenet.

3. **Potential Evaluation (Modern Fortran)**

- Uses trained models to compute potential energy and atomic forces.

- Designed for high-performance molecular dynamics and structure optimization.
- 
**Features**

- **Modern Fortran implementation** for descriptor generation and potential evaluation

- **Training compatibility with aenet** (until the Julia module is released)


**Installation**

**Prerequisites**

• A modern Fortran compiler (GFortran, Intel Fortran, or NVIDIA Fortran)

• CMake (version 3.15 or later)


**Build Instructions**

```
git clone https://github.com/cometscome/BPNET.git
cd BFNET
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**License**

BFNET is released under the MIT License. See the [LICENSE](LICENSE) file for details.


**Contributors**

- [Yuki Nagai](https://github.com/cometscome)

- Contributions welcome! Feel free to submit issues and pull requests.

  
**Acknowledgments**


BFNET is inspired by [aenet](http://ann.atomistic.net) and aims to provide a modern, maintainable alternative.
