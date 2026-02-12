# Replication of "Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model"
A revised version of the code used for the simulation study in "Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model" by Dr. [Karsten Schweikert](https://github.com/karstenschweikert) and I.

## Introduction

This repository contains the replication code for the following working paper:

Franjic, Domenic and Schweikert, Karsten, Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model (October 30, 2024), Last revised: 12 Feb 2026. Available at SSRN: [https://ssrn.com/abstract=4733872](https://ssrn.com/abstract=4733872) or [https://ssrn.com/abstract=4733872](http://dx.doi.org/10.2139/ssrn.4733872) 

## SimulationStudyReplication.cpp

This file, together with the routines in ``./Internals/``, provides the code to replicate the simulation results of our study.

### Features

- **Fast and Parallelised Cross-Validation**: Implements a parallelised random hyper-parameter search for efficient cross-validation.
- **Flexible Simulation Parameters**: Allows a high degree of customisation of the model parameterisation.
- **Compatibility**: Works on modern Windows and Linux operating systems.
- **Open-Source**: Distributed under the GNU General Public License v3.0.

### Prerequisites

- **Eigen** (version 3.4.0 or later): A C++ template library for linear algebra. [Eigen Website](https://eigen.tuxfamily.org/)
- **OpenMP-compatible C++ Compiler**: Such as GCC (version 5.0 or later) or MSVC (Visual Studio 2019 or later).

### Installation

1. Make sure that a compatible version of Eigen3 (3.4.0 or later) is installed.
2. Install a OpenMP compatible C++ compiler such as GNU C++ Compiler or MSVC.
3. Clone the Repository
   ```bash
   git clone https://github.com/yourusername/ReplicationNowcastingMacroVarsWithSDFM.git
   cd ReplicationNowcastingMacroVarsWithSDFM
4. The repository is ready to compile. Using MinGW-w64 GCC with reasonable speed optimisation on windows, a potential compiler call could look something like this:
   ```bash
   g++ -std=c++14 -Wall -O3 -march=native -mfma -fopenmp -DNDEBUG -I "C:\Path\To\Eigen" SimulationStudyReplication.cpp Internals\*.cpp -o SimulationStudyReplication.exe
   ```
   On Linux, a possible compiler call might be
   ```bash
   g++ -std=c++14 -Wall -O3 -march=native -mfma -fopenmp -DNDEBUG -I /path/to/eigen3 SimulationStudyReplication.cpp Internals/*.cpp -o SimulationStudyReplication
   ```
   Note that the speed optimisations are highly encouraged due to the computational complexity of the cross-validation scheme.

### Usage

Most of the model parameters are set at compile time. However, it is possible to interrupt the simulations and restart them at a later time. To re-parameterise the simulation study, it is generally sufficient to change the "hard" parameters at the beginning of SimulationStudyReplication.cpp. Further or more general changes to the model parameterisation require a reformulation of the structure in the data-generating part of SimulationStudyReplication.cpp or at deeper levels.

## EmpiricalStudyReplication.r

This file provides the code to replicate the empirical results of our study.

### Prerequisites

- **TwoStepSDFM**: A beta version of my ``R`` package that implements, among other things, the estimation, cross-validation, and nowcasting schemes outlined in our study. [TwoStepSDFM Github Repo](https://github.com/SiSanchopancho/TwoStepSDFM.git)

### Usage

Before using the code, you must download all available FRED-MD vintages from the  [FRED webiste](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) (McCracken, M. W. 2024. “FRED-MD and FRED-QD: Monthly and Quarterly Databases for Macroeconomic Research.” Federal Reserve Bank of St. Louis). You also need a single FRED-QD file to extract the transformation code corresponding to the US GDP level series. We recommend using the ``fred-qd_2024m12.csv`` dataset, since it is already referenced in the code. 

Using the script is straightforward. The first time you run it, uncomment and execute the data‐download block to load, preprocess, and restructure the mixed frequency vintage files. After that step completes, you can execute the script without further modification.

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

© 2024-2026 Domenic Franjic

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This work is partially based on the LARS-EN and SPCA algorithms found in:

- Zou, H., Hastie, T., & Tibshirani, R. (2006). *Sparse Principal Component Analysis*. Journal of Computational and Graphical Statistics, 15(2), 265-286.
- Normal Splines. (2019, February 25). Algorithms for updating the Cholesky factorization. Normal Splines Blog. (https://normalsplines.blogspot.com/2019/02/algorithms-for-updating-cholesky.html)

I also utilise the following libraries:

- **Eigen 3**: Guennebaud, G., Jacob, B., & Others (2010). *Eigen v3*. [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org)
- **OpenMP**: OpenMP Architecture Review Board (2020). *OpenMP Application Programming Interface Version 5.1*. [https://www.openmp.org/spec-html/5.1/](https://www.openmp.org/spec-html/5.1/)

## Contributing

As this repo provides the code to replicate the results of our study, it is not possible to contribute.

## Support

If you have any questions or need assistance, please open an issue on the GitHub repository or contact us via email.

## Contact

- **Name**: Domenic Franjic
- **Institution**: University of Hohenheim
- **Department**: Econometrics and Statistics, Core Facility Hohenheim
- **E-Mail**: franjic@uni-hohenheim.de
