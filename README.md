
# **BPNET: A Modern Fortran Machine Learning Potential Code**

BPNET (**Behler-Parrinello neural NETwork potential**) is an open-source, modern Fortran implementation of machine learning potential construction. It is inspired by [aenet](http://ann.atomistic.net) and aims to provide a more maintainable and extensible alternative while introducing new features and optimizations.


**Note:** BPNET is still under development and does not yet support all features of aenet. However, it is designed to **extend beyond aenet** by incorporating new capabilities.

**Architecture**

BPNET is structured into **three key components**, leveraging **modern Fortran** for high-performance computation and **Julia** (under development) for flexible model training:

1. **Descriptor Generation (Modern Fortran)**

- Computes atomic environment descriptors based on Behler-Parrinello-type symmetry functions. Now only Chebyshev descriptor (compatible with aenet v2.0.4 not 2.0.3) is supported.

2. **Neural Network Training (Julia, in development)**

- The training module is currently under development and will be released in a separate repository.

- In the meantime, **aenet can be used for training**, as BPNET’s descriptor format is compatible with aenet.

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
cd BPNET
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Descriptor Generation**
BPNET generates atomic environment descriptors using Behler-Parrinello-type symmetry functions. The current implementation supports Chebyshev descriptors, which are compatible with aenet v2.0.4 and v2.0.3.

For aenet v2.0.4, the descriptor generation is compatible with the following settings:
```
BASIS type=Chebyshev
radial_Rc = 8.0  radial_N = 16 angular_Rc = 6.5  angular_N = 4 version=1
```

For aenet v2.0.3, the descriptor generation is compatible with the following settings:
```
BASIS type=Chebyshev
radial_Rc = 8.0  radial_N = 16 angular_Rc = 6.5  angular_N = 4 
```

**License**

BPNET is released under the MIT License. See the [LICENSE](LICENSE) file for details.


**Contributors**

- [Yuki Nagai](https://github.com/cometscome)

- Contributions welcome! Feel free to submit issues and pull requests.

  
**Acknowledgments**


BPNET is inspired by [aenet](http://ann.atomistic.net) and aims to provide a modern, maintainable alternative.
