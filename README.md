# ReplicationNowcastingMacroVarsWithSDFM
A revised version of the code used for the simulation study in "Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model" by Dr. [Karsten Schweikert](https://github.com/karstenschweikert) and I.

## Introduction

This repository contains the replication code for the following working paper:

Franjic, Domenic and Schweikert, Karsten, *Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model* (February 21, 2024). Available at SSRN: https://ssrn.com/abstract=4733872 or http://dx.doi.org/10.2139/ssrn.4733872

In short, the code draws data from a dynamic factor model, i.e., a linear Gaussian state-space model. Different degrees of measurement-error cross-correlation, degrees of sparsity in the loadings matrix, and different numbers of factors are considered. Each parametrisation is repeated 1000 times. The sparse estimators are validated using random hyper-parameter search in combination with either the BIC or time-series cross-validation. The validation schemes are parallelised to achieve a reasonable computational performance for larger parametrisations. Finally, the code will provide the *mean squared nowcasting error* (MSNE). Fore more details, see our paper.

## Features

- Provides a fast and parallelised implementation for the cross-validation scheme of the sparse two-step SDFM estimator.
- The program is compatible with all modern Windows and Linux OS.

## Prerequisites

- [Eigen](https://eigen.tuxfamily.org/) (version 3.4.0 or later)
- OpenMP compatible C++ compiler.

## Installation

1. Make sure that a compatible version of Eigen3 (3.4.0 or later) is installed.
2. Install a OpenMP compatible C++ compiler such as GNU C++ Compiler or MSVC.
3. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/yourproject.git
4. The repository is ready to compile. Using gcc with reasonable speed optimisation on windows, a potential compiler call could look something like this:
   ```bash
   @echo off
   g++ -std=c++14 -Wall -g -O3 -march=native -mfma -fopenmp -DNDEBUG -IPath\To\Eigen" ReplicationNowcastingMacroVarsWithSDFM.cpp Internals\Filtering.cpp Internals\DataGen.cpp Internals\DataHandle.cpp Internals\SparsePCA.cpp Internals\CholUpDown.cpp Internals\SparseDFM.cpp Internals\CrossVal.cpp Internals\Forecast.cpp -o ReplicationNowcastingMacroVarsWithSDFM
   ```

## License

© 2024 Domenic Franjic

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

This work is partially based on the LARS-EN and SPCA algorithm found in:

- Zou, H., Hastie, T., & Zou, M. H. (2016). *Package "elasticnet"*. Available at https://cran.r-project.org/web/packages/elasticnet/index.html

This work makes use of the C++ libraries Eigen3 and OpenMP:

- Guennebaud, Ga\"{e}l,  Jacob, Beno\^{i}t et al. (2010). *Eigen v3*. Available at http://eigen.tuxfamily.org
- OpenMP Architecture Review Board. (2020). *OpenMP Application Programming Interface Version 5.1.*. Available at https://www.openmp.org/wp-content/uploads/OpenMP-API-Specification-5-1.pdf
