# ReplicationNowcastingMacroVarsWithSDFM
A revised version of the code used for the simulation study in "Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model" by Dr. [Karsten Schweikert](https://github.com/karstenschweikert) and I.

## Introduction

This repository contains the replication code for the following working paper:

Franjic, Domenic and Schweikert, Karsten, *Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model* (February 21, 2024). Available at SSRN: https://ssrn.com/abstract=4733872 or http://dx.doi.org/10.2139/ssrn.4733872

In short, the code draws data from a dynamic factor model, i.e., a linear Gaussian state-space model, according to
$$
\mathbf{x}_{t}^* &= \boldymbol{\Lambda}\mathbf{f}_t + \bm{\xi}_t^*
\bm{\Phi}^P(\mathbb{L})\mathbf{f}_{t} &= \bm{\epsilon}_t
\mathbb{E}(\bm{\xi}_t^*{\bm{\xi}_t^*}') &= (\sigma_{n,m})_{n=1,m=1}^{N,N}=:\bm{\Sigma},
$$
for $t=1/3,2/3,\dots,T$, where $\bm{\epsilon}_t\sim\mathcal{N}(\mathbf{0}_R,\bm{\Omega})$ represents the state error

## Features

- Feature 1
- Feature 2

## Prerequisites

- [Eigen](https://eigen.tuxfamily.org/) (version X.X or later)
- [OpenBLAS](https://www.openblas.net/) (version X.X or later)

## Installation

1. Install Eigen and OpenBLAS.
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/yourproject.git

## License

© 2024 Domenic Franjic

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

This work is partially based on thefollowing work

- Zou, H., Hastie, T., & Zou, M. H. (2016). *Package "elasticnet"*. Available on https://cran.r-project.org/web/packages/elasticnet/index.html

This work makes use of the following C++ linear algebra libraries

- Guennebaud, Ga\"{e}l,  Jacob, Beno\^{i}t et al. (2010). *Eigen v3*. Available at http://eigen.tuxfamily.org
