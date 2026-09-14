# Replication of "Nowcasting Macroeconomic Variables with a Sparse Mixed-Frequency Dynamic Factor Model"
A revised version of the code used for the simulation study in "Nowcasting Macroeconomic Variables with a Sparse-Mixed Frequency Dynamic Factor Model" by Dr. [Karsten Schweikert](https://github.com/karstenschweikert) and I.

## Introduction

This repository contains the replication code for the following working paper:

Franjic, Domenic and Schweikert, Karsten, Nowcasting Macroeconomic Variables with a Sparse Mixed-Frequency Dynamic Factor Model (October 30, 2024), Last revised: 13 July 2026. Available at SSRN: [https://ssrn.com/abstract=4733872](https://ssrn.com/abstract=4733872) or [https://ssrn.com/abstract=4733872](http://dx.doi.org/10.2139/ssrn.4733872) 

## Repository structure

### Structure

```text
ReplicationNowcastingMacroVarsWithSDFM/
├── Internals/                         # C++ source and header files used by the simulation study
├── RHelper/                           # R helper functions for evaluating the empirical results
├── TwoStepSDFM_0.3.0.3.tar.gz         # Exact package version used to generate the empirical results
├── EmpiricalStudyReplication.R        # Main script for reproducing the empirical results including all necessary about data extraction
├── SimulationStudyReplication.cpp     # Main script for reproducing the simulation results 
├── README.md
├── LICENSE
```

### Code and manuscript outputs

| Manuscript output | Script | Notes |
|---|---|---|
| Table 1 | `SimulationStudyReplication.cpp` | The simulation study requires approximately 2–3 weeks to run, depending on the computing architecture. |
| Figure 1 & Table 2 | `EmpiricalStudyReplication.R` | Reproduces the empirical results. |
| Figures 2 & 3 | `GroupImportance.R` | Produces the group-importance figures. |

## SimulationStudyReplication.cpp

This file, together with the routines in ``./Internals/``, provides the code to replicate the simulation results of our study.

### Prerequisites

#### C++ computing environment

- **Platform:** x86_64
- **Operating system:** Ubuntu 22.04
- **Linux kernel:** 5.15.0-190-generic
- **CPU:** Intel Xeon Gold 6348H @ 2.30 GHz
- **Logical CPUs:** 24
- **Compiler:** g++ 11.4.0
- **Parallelisation:** OpenMP with 24 threads

#### C++ library dependencies

| Library | Version |
|---|---|
| Eigen | 3.4.0 |

**Note:** An OpenMP-compatible C++ Compiler, Such as GCC (version 5.0 or later) or MSVC (Visual Studio 2019 or later) is required.

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

#### R computing environment

- **R version:** 4.3.1 (2023-06-16)
- **Platform:** x86_64-w64-mingw32/x64 (64-bit)
- **Operating system:** Windows 11 x64 (build 26200)

#### R package dependencies

| Package | Version |
|---|---|
| rstudioapi | 0.18.0 |
| zoo | 1.8-15 |
| lubridate | 1.9.5 |
| BVAR | 1.0.5 |
| readxl | 1.5.0 |
| alfred | 0.2.1 |
| TwoStepSDFM | 0.3.0.3 |
| sandwich | 3.1-1 |
| lmtest | 0.9-40 |
| tidyr | 1.3.2 |
| dplyr | 1.2.1 |
| car | 3.1-5 |
| stargazer | 5.2.3 |
| murphydiagram | 0.12.2 |
| dynlm | 0.3-6 |
| modelsummary | 2.6.0 |
| ggplot2 | 4.0.3 |
| xtable | 1.8-8 |

**Note on TwoStepSDFM**: To replicate the results, version 0.3.0.3 of my `R` package, which implements, among other things, the estimation, cross-validation, and nowcasting schemes outlined in our study, is required. As of now, this version of the package is not available on CRAN. A .tar.gz-ball of version 0.3.0.3 is found in this repo. For the most recent version of the package, see the [TwoStepSDFM GitHub repository](https://github.com/SiSanchopancho/TwoStepSDFM.git). For the most stable version of the package see [TwoStepSDFM on CRAN](https://cran.r-project.org/web/packages/TwoStepSDFM/index.html). Please note that using the current version of the package may produce results that differ slightly from those presented in the study. 

### Data

Before using the code, you must download all available FRED-MD vintages from the  [FRED webiste](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) (McCracken, M. W. 2024. “FRED-MD and FRED-QD: Monthly and Quarterly Databases for Macroeconomic Research.” Federal Reserve Bank of St. Louis). You also need a single FRED-QD file to extract the transformation code corresponding to the US GDP level series. We recommend using the ``fred-qd_2024m12.csv`` dataset, since it is already referenced in the code. 

Using the script is straightforward. The first time you run it, uncomment and execute the data‐download block to load, preprocess, and restructure the mixed frequency vintage files. After that step completes, you can execute the script without further modification.

## Expected runtime

- **Simulation study (Table 1):** approximately 2–3 weeks, depending on the computing architecture.
- **Empirical analysis (Figure 1 and Table 2):** 5-6 hourse, depending on the computing architecture.
- **Group-importance analysis (Figures 2–3):** ~1 hour, depending on the computing architecture.

## Simulation reproducibility

The simulation code uses a fixed random seed (`18092024`). The simulations were compiled using CPU-specific optimisation (`-march=native`), fused multiply-add instructions (`-mfma`), and OpenMP parallelisation. Consequently, reruns on different hardware or software environments may yield small numerical deviations and are not guaranteed to be bitwise identical.

# Output files

The repository contains the complete source code required to generate all manuscript outputs. Output files are generated by running the scripts listed in the “Code and manuscript outputs” section. Pre-generated output files are not included in the repository.

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
- **Reproducibility package assembled:** 14. September 2026
