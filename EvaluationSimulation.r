# SPDX-License-Identifier: GPL-3.0-or-later #
#
# Copyright (C) 2024-2026 Domenic Franjic
#
# This file is part of ReplicationNowcastingMacroVarsWithSDFM.
#
# ReplicationNowcastingMacroVarsWithSDFM is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ReplicationNowcastingMacroVarsWithSDFM is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ReplicationNowcastingMacroVarsWithSDFM. If not, see <https://www.gnu.org/licenses/>.

library(rstudioapi)
library(ggplot2) 
library(reshape2) 
library(xtable)
library(ggpubr)
library(tidyr)

rm(list = ls())

setwd(dirname(getActiveDocumentContext()$path))

# Set up parameters of the current simulation set-up aimed to be investigated
no_of_observations <- 303
no_of_diff_factors <- 3
no_of_models <- 3
diff_corr_simulated <- c(4212, 1737, 995, 500)
avg_abs_corr_simulated <- c(0.2, 0.3, 0.4, 0.5)
no_of_vars_checked <- c(25, 50)
no_of_factors_checked <- c(1, 2, 3)
no_of_diff_corr <- length(diff_corr_simulated)
diff_zero_prop_simulated <- c(0, 40, 80)
no_of_diff_prop <- length(diff_zero_prop_simulated)

msne <- matrix(NaN, no_of_diff_prop * no_of_diff_corr * length(no_of_vars_checked) * length(no_of_factors_checked),
               8)
colnames(msne) <- c("D Non-Diag", "BIC Non-Diag", "CV Non-Diag", 
                    "D Diag", "BIC Diag", "CV Diag", 
                    "ARMA", "UM")
rownames(msne) <- rep("DUMMY", no_of_diff_prop * no_of_diff_corr * length(no_of_vars_checked) * length(no_of_factors_checked))

# Start loop over the number of variables #

row_ind <- 1
for(no_of_variables in no_of_vars_checked){
  
  # Start loop over number of factors #
  
  parametrisation_ind <- 1
  plot_list <- list()
  for (no_of_factors in no_of_factors_checked) {
    
    # Create temporal storage matrices for storing the data
    mean_reduction <- matrix(NaN, no_of_diff_prop, no_of_diff_corr * no_of_models)
    colnames(mean_reduction) <- rep(paste0("beta = ", round(diff_corr_simulated / 1000, 2)), no_of_models)
    rownames(mean_reduction) <- paste0("p = ", round(diff_zero_prop_simulated / 100, 2))
    
    success_prob <- matrix(NaN, no_of_diff_prop, no_of_diff_corr * no_of_models)
    colnames(success_prob) <- rep(paste0("beta = ", round(diff_corr_simulated / 1000, 2)), no_of_models)
    rownames(success_prob) <- paste0("p = ", round(diff_zero_prop_simulated / 100, 2))
    
    # Start loop over the different proportions of zeros in the loading matrix #
    
    for (curr_prop in diff_zero_prop_simulated) {
      
      # Start loop over the different degrees of correlaiton of the measurement error matrix #
      
      for (curr_corr in diff_corr_simulated) {
        
        # Load the data
        name_append <- paste0(
          "_K", no_of_factors, "_N", no_of_variables, "_T", no_of_observations, "_p", curr_prop, "_b", curr_corr, "_c", 1, "_H", 3, "_TCV", no_of_observations - 3, "_G1", 25 * no_of_factors,
          "_seed", 18092024
        )
        data_loaded_check <- tryCatch(
          {
            MSFE <- read.csv(paste0("./SimulationResults/diag_non_diag_uniform_lambda_results_T_", no_of_observations, "_N_", no_of_variables, "_K_", no_of_factors, "/RMSFE", name_append, ".csv"), header = FALSE)
            FALSE
          },
          error = function(e) {
            warning(paste0("dataset: ", name_append, " not loaded."))
            return(TRUE)
          }
        )
        if(data_loaded_check){
          next
        }
        
        msne[row_ind, ] <- colMeans(MSFE, na.rm = TRUE)
        rownames(msne)[row_ind] <- paste0("vars_", no_of_variables, "_fact_", no_of_factors, "_prop_", curr_prop, "_corr_", curr_corr)
        row_ind <- row_ind + 1
      }
      
      # Start loop over the different degrees of correlaiton of the measurement error matrix #
      
    }
    
    # End loop over the different proportions of zeros in the loading matrix #
    
  }
  
  # End loop over number of factors #
  
}

# End loop over the number of variables #

# Create result tables #

makeComparinsTable <- function(model, benchmark){
  rel_performance <- as.data.frame(1 - model / benchmark)
  rel_performance$dimensions <- sub("^(.*fact_[0-9]+)_(.*)$", "\\1", rownames(rel_performance))
  rel_performance$paramet <- sub("^(.*fact_[0-9]+)_(.*)$", "\\2", rownames(rel_performance))
  results_formatted <- rel_performance[1:12, 3, drop = FALSE]
  running_ind <- 1
  for(ind in seq(1, dim(rel_performance)[1] - 8, 12)){
    results_formatted <- merge(results_formatted, rel_performance[ind:(ind + 11), c(1, 3), drop = FALSE], by = "paramet")
    colnames(results_formatted)[running_ind + 1] <- rel_performance[ind, 2]
    running_ind <- running_ind + 1
  }
  rownames(results_formatted) <- results_formatted$paramet
  results_formatted <- results_formatted[, 2:7]
  results_formatted
}

# Create relative comparison tables for CV Non-Diag against D Non-Diad
results_formatted <- makeComparinsTable(msne[, 3, drop = FALSE], msne[, 1, drop = FALSE])
results_formatted

# Create relative comparison tables for BIC Non-Diag against D Non-Diad
results_formatted_bic <- makeComparinsTable(msne[, 2, drop = FALSE], msne[, 1, drop = FALSE])
results_formatted_bic

# Create relative comparison tables for CV Non-Diag against D Diad
results_formatted_ndiag_diag <- makeComparinsTable(msne[, 3, drop = FALSE], msne[, 5, drop = FALSE])
results_formatted_ndiag_diag

# Create relative comparison tables for BIC Non-Diag against D Diad
bic_results_formatted_ndiag_diag <- makeComparinsTable(msne[, 2, drop = FALSE], msne[, 5, drop = FALSE])
bic_results_formatted_ndiag_diag

# Make tex tables #

cellColourValues <- function(value){
  paste0("\\g{", sprintf("%.2f", 100 * value), "}")
}


makeTexTable <- function(label, table){
  cat(
    "% ", as.character(Sys.time()), "
    \\begin{table}[ht]
    \\begin{center}\n",
    label,
    "\\resizebox{\\textwidth}{!}{\\begin{tabular}{rrcccccc}
    \\hline
    \\hline
     & & \\multicolumn{3}{c}{$N=25$} & \\multicolumn{3}{c}{$N = 50$}\\\\
     & & $R = 1$ & $R = 2$ & $R = 3$ & $R = 1$ & $R = 2$ & $R = 3$ \\\\ 
    \\hline
    \\hline\n")
  
  # Results for s = 0
  cat(
    "\\multirow{4}{*}{$s = 0.0$} & $\\bar{\\rho}\\approx 0.2$ & ",
    paste0(apply(table[1, ], 1, cellColourValues), collapse = " & "),
    "\\\\\n"
  )
  correlation_vector <- c(0.3, 0.4, 0.5)
  ind <- 1
  for(row in 2:4){
    cat(
      paste0(" & $\\bar{\\rho}\\approx", correlation_vector[ind], "$ & "),
      paste0(apply(table[row, ], 1, cellColourValues), collapse = " & "),
      "\\\\\n"
    )
    ind <- 1 + ind
  }
  
  # Results for s = 0.4
  cat(
    "\\hline
    \\multirow{4}{*}{$s = 0.4$} & $\\bar{\\rho}\\approx 0.2$ & ",
    paste0(apply(table[5, ], 1, cellColourValues), collapse = " & "),
    "\\\\\n"
  )
  ind <- 1
  for(row in 6:8){
    cat(
      paste0(" & $\\bar{\\rho}\\approx", correlation_vector[ind], "$ & "),
      paste0(apply(table[row, ], 1, cellColourValues), collapse = " & "),
      "\\\\\n"
    )
    ind <- 1 + ind
  }
  
  # Results for s = 0.8
  cat(
    "\\hline
    \\multirow{4}{*}{$s = 0.8$} & $\\bar{\\rho}\\approx 0.2$ & ",
    paste0(apply(table[9, ], 1, cellColourValues), collapse = " & "),
    "\\\\\n"
  )
  ind <- 1
  for(row in 10:12){
    cat(
      paste0(" & $\\bar{\\rho}\\approx", correlation_vector[ind], "$ & "),
      paste0(apply(table[row, ], 1, cellColourValues), collapse = " & "),
      "\\\\\n"
    )
    ind <- 1 + ind
  }

  # Table bottom
  cat("\\hline
  \\hline
  \\end{tabular}}
  \\end{center}
  \\begin{tablenotes}
  \\scriptsize\\item Note: $\\rho$ refers to the average absolute correlation-coefficients of the measurement errors. $s$ refers to the degree of sparsity which indicates the proportion of zero elements in each column of the factor loading matrix $\\bm{\\Lambda}$. The colour coding indicates the relative performance of the sparse models versus the benchmark per number of factors, with green indicating better and red indicating worse performance. The gradient of the colour indicates the magnitude of the MSNE improvement/deterioration.
  \\end{tablenotes} 
  \\end{table}")
}

# Create relative comparison tables for CV Non-Diag against D Non-Diad
makeTexTable("\\caption{MSNE reduction in \\% over 1000 simulated nowcasting exercises using $T = 100$ observations with Hyper-Parameters validated via CV using a non-diagonal variance-covariance matrix for the dense and sparse two-step estimator\\label{tbl::T100_msne_reduction_non_diag}}", 
             results_formatted)

# Create relative comparison tables for BIC Non-Diag against D Non-Diad
makeTexTable("\\caption{MSNE reduction in \\% over 1000 simulated nowcasting exercises using $T = 100$ observations with Hyper-Parameters validated via the BIC using a non-diagonal variance-covariance matrix for the dense and sparse two-step estimator\\label{tbl::T100_msne_reduction_non_diag_bic}}", 
             results_formatted_bic)

# Create relative comparison tables for BIC Non-Diag against D Diad
makeTexTable("\\caption{MSNE reduction in \\% over 1000 simulated nowcasting exercises using $T = 100$ observations with Hyper-Parameters validated via the BIC using a non-diagonal variance-covariance matrix for the sparse and a diagonal variance covariance matrix for the dense two-step estimator\\label{tbl::T100_msne_reduction_non_diag_diag}}", 
             bic_results_formatted_ndiag_diag)

# Create relative comparison tables for CV Non-Diag against D Diad
makeTexTable("\\caption{MSNE reduction in \\% over 1000 simulated nowcasting exercises using $T=100$ observations\\label{tbl::T100_msne_reduction}}", 
             results_formatted_ndiag_diag)

