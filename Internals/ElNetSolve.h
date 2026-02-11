/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2026 Domenic Franjic
 *
 * This file is part of ReplicationNowcastingMacroVarsWithSDFM.
 *
 * ReplicationNowcastingMacroVarsWithSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * ReplicationNowcastingMacroVarsWithSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ReplicationNowcastingMacroVarsWithSDFM. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef ENET_SOLVER
#define ENET_SOLVER

 // Including external libraries
#include <Eigen>
#include <stdlib.h>
#include <math.h>
#include <cfloat>

// Including internal libraries
#include "DataHandle.h"
#include "CholUpDown.h"

/* LARS Algorithmus */
/* Sources:
- Zou, H., Hastie, T., & Zou, M. H. (2016). Package ‘elasticnet’. (https://cran.r-project.org/web/packages/elasticnet/index.html)
- Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least angle regression. (DOI: 10.1214/009053604000000067)
- Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society Series B: Statistical Methodology, 67(2), 301-320. (DOI: https://doi.org/10.1111/j.1467-9868.2005.00527.x)
*/
template<bool return_penalties>
Eigen::MatrixXd LARS(
  const Eigen::VectorXd& y,  // VOI
  const Eigen::MatrixXd& X, // Predictors
  const double& l2 = 10e-6, // l2 Penalty
  double l1 = NAN, // l1 penalty (used for stopping)
  int selected_in = INT_MAX, // Number of selected variables
  int steps = INT_MIN, // Number of steps until stopping
  const double& comp_null = 10e-15// Computational zero
)
{

  /* Dummies */

  // Integers
  const int N = static_cast<int>(X.cols()), T = static_cast<Eigen::Index>(X.rows());
  int selected = selected_in, just_left = -1, curr_var_index = -1, gamma_index = -1, size_A_C = N, size_A = 0, s = 0, drop_index = -1;

  // Bools
  bool drop = 0;

  // Reals
  double gamma = DBL_MAX, gamma_hat = DBL_MAX, A_A = DBL_MAX, gamma_tilde = DBL_MAX, C_hat = DBL_MAX, l2_sqrt = sqrt(l2), l2_sqrt_inv = 1 / sqrt(1 + l2);

  // Matrices
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(N, N);

  // Vectors
  Eigen::VectorXd w_A_vec = Eigen::VectorXd::Zero(N), L_T_G_inv = Eigen::VectorXd::Zero(N), G_A_inv_one = Eigen::VectorXd::Zero(N), g_tilde = Eigen::VectorXd::Zero(N),
    beta_curr = Eigen::VectorXd::Zero(N), u_A = Eigen::VectorXd::Zero(N + T), a(N), residuals = Eigen::VectorXd::Zero(T + N),
    c = l2_sqrt_inv * X.transpose() * y, beta_hat = Eigen::VectorXd::Zero(N), sign = Eigen::VectorXd::Zero(N);

  Eigen::VectorXi A_C_set = Eigen::VectorXi::LinSpaced(N, 0, N - 1), A_set = Eigen::VectorXi::Constant(N, -1);

  // Special cases and misshandling
  if (!(std::isnan(l1)) && (selected != INT_MAX))
  {
    std::cout << '\n' << "Error! Both penalty and number of variables is provided. Termination criterion is ambiguous. Only provide either of them." << '\n';
    return EXIT_FAILURE * Eigen::VectorXd::Ones(N);
  }
  else if (l1 < 0 || (selected < 0 && selected != INT_MAX) || (steps < 0 && steps != INT_MIN))
  {
    std::cout << "l1 = " << l1 << '\t' << "selected = " << selected << '\t' << "steps = " << steps << '\t';
    std::cout << '\n' << "Error! A termination criterion ('l1', 'selected', or 'steps') has negative value." << '\n';
    return EXIT_FAILURE * Eigen::VectorXd::Ones(N);
  }
  else if ((!(std::isnan(l1)) && ((2 * c.cwiseAbs().col(0).maxCoeff() / l2_sqrt_inv) <= l1)) || selected == 0)
  {
    return Eigen::VectorXd::Zero(N);
  }

  // Set the default stopping criterion with respect to number of variables selected
  if (l2 == 0.)
  {
    selected = (N <= (T - 1)) ? N : T - 1;
  }
  else if (N < selected)
  {
    selected = N;
  }

  // Redefinitions
  steps = ((steps == INT_MIN) ? 50 * (T <= N - 1 ? T : (N - 1)) : steps);
  residuals.head(T) = y;
  residuals.tail(N).setZero();
  Eigen::VectorXd pen(steps + 1);
  pen(0) = c.cwiseAbs().col(0).maxCoeff();

  /* LARS-EN loop */
  while (s < steps && size_A < selected)
  {

    // Calculate current maximum correlation and retrieve the index of the corresponding variable evaluated at the current step
    C_hat = c(A_C_set.head(size_A_C)).cwiseAbs().maxCoeff(&just_left);

    if (drop == 0)
    {

      /* Case 1: Previously no variable has just been dropped */

      // Add the variable to the active set and rase the variable from the inactive set
      curr_var_index = A_C_set(just_left);
      --size_A_C;
      ++size_A;
      A_set(size_A - 1) = curr_var_index;
      sign(size_A - 1) = double(0 < c(curr_var_index)) - double(c(curr_var_index) < 0);
      A_C_set(seq(just_left, N - 2)) = A_C_set(seq(just_left + 1, N - 1)).eval();
      A_C_set(N - 1) = -1;

      if (s == 0)
      {

        // Create Gramm-Matrix
        L(0, 0) = sqrt(((X.col(curr_var_index).transpose() * X.col(curr_var_index))(0, 0) + l2) / (1 + l2));

      }
      else
      {

        // Update Cholesky
        CholUpDown::cholUpdate(L, X(all, curr_var_index), X(all, A_set.head(size_A - 1)), -1, l2);

      }
    }
    else if (drop == 1 || size_A == N)
    {

      /* Case 1: Previously no variable has just been dropped */

      // No new variable will be added since the last step has led to taking an "incomplete" LARS step
      curr_var_index = -1;
      drop = 0;

    }

    // Calculate equiengular vector
    G_A_inv_one.head(size_A) = sign.head(size_A);
    L.topLeftCorner(size_A, size_A).triangularView<Eigen::Lower>().solveInPlace(G_A_inv_one.head(size_A));
    L.topLeftCorner(size_A, size_A).transpose().triangularView<Eigen::Upper>().solveInPlace(G_A_inv_one.head(size_A));
    A_A = 1.0 / (std::sqrt((G_A_inv_one.head(size_A).transpose() * sign.head(size_A))(0)));
    Eigen::VectorXd w_A = A_A * G_A_inv_one.head(size_A);
    w_A_vec(A_set.head(size_A)) = w_A;
    u_A.head(T) = l2_sqrt_inv * (X(all, A_set.head(size_A)) * w_A);
    u_A.tail(N) = l2_sqrt * l2_sqrt_inv * w_A_vec;

    /* Computing the step-sizes */

    // Calculate the maximum feasible step size (LARS-LASSO-Modification)
    gamma_tilde = DBL_MAX;
    for (int nn : A_set.head(size_A))
    {

      double g_curr = (-1 * (beta_hat(nn) / w_A_vec(nn)));
      if (g_curr < gamma_tilde && comp_null < g_curr)
      {

        gamma_tilde = g_curr;
        gamma_index = nn;

      }

    }

    if (N == size_A)
    {
      // For the last step, just go all the way and set gamma_hat to the maximum correlation
      gamma_hat = C_hat / A_A;

    }
    else
    {

      // Compute step size of the current LARS step
      a.head(size_A_C) = (X(all, A_C_set.head(size_A_C)).transpose() * u_A.head(T) + l2_sqrt * u_A.tail(N)(A_C_set.head(size_A_C))) * l2_sqrt_inv;
      gamma_hat = C_hat / A_A;

      for (int nn = 0; nn < size_A_C; ++nn)
      {

        double CAm = (C_hat - c(A_C_set(nn))) / (A_A - a(nn));
        double CAp = (C_hat + c(A_C_set(nn))) / (A_A + a(nn));

        if (CAm < gamma_hat && comp_null < CAm)
        {
          gamma_hat = CAm;
        }
        else if (CAp < gamma_hat && comp_null < CAp)
        {
          gamma_hat = CAp;
        }
      }
    }

    /* Updating beta_hat, the residuals, and the correlation vector */

    // If gamma_hat is bigger then gamma_tilde, it is not possible to do a complete step due to sign restrictions on the coefficients
    // In this case the variable that would switch signs is dropped and the step is only as large as necessary to make the coefficient
    // of the variable that is dropped go to zero.
    gamma = (gamma_tilde < gamma_hat) ? gamma_tilde : gamma_hat;
    beta_curr = beta_hat;
    beta_hat += (gamma * w_A_vec);
    residuals -= (gamma * u_A);
    c = (X.transpose() * residuals.head(T) + l2_sqrt * residuals.tail(N)) * l2_sqrt_inv;

    // Check whether the program should be stopped early due to the li penalty
    pen(s + 1) = pen(s) - std::abs(gamma * A_A);
    if (!(std::isnan(l1)) && (((pen(s + 1) * 2) / l2_sqrt_inv) <= l1)) {

      double ps1_l2 = pen(s + 1) * 2 / l2_sqrt_inv;
      double 	ps_l2 = pen(s) * 2 / l2_sqrt_inv;
      beta_hat = (((ps_l2 - l1) / (ps_l2 - ps1_l2)) * beta_hat + ((l1 - ps1_l2) / (ps_l2 - ps1_l2)) * beta_curr) * l2_sqrt_inv;

      // Return the penalties or the coefficient vector depending on the needs
      if (return_penalties) {
        return pen(seq(0, s)) * (2.0 / l2_sqrt_inv);
      }
      else {
        return beta_hat;
      }

    }

    /* Add or remove variacles */

    if (gamma_tilde < gamma_hat)
    {

      /* Drop situation (for computational reasons cast the coefficients as zeros explicitly) */

      drop = 1;
      beta_hat(gamma_index) = 0.;
      w_A_vec(gamma_index) = 0.;
      u_A.tail(N)(gamma_index) = 0.;
      ++size_A_C;
      A_C_set(size_A_C - 1) = gamma_index;
      (A_set.array() == gamma_index).maxCoeff(&drop_index);
      if (drop_index != N - 1)
      {

        A_set(seq(drop_index, N - 2)) = A_set(seq(drop_index + 1, N - 1)).eval();
        sign(seq(drop_index, N - 2)) = sign(seq(drop_index + 1, N - 1)).eval();

      }
      A_set(N - 1) = -1;
      sign(N - 1) = 0;
      --size_A;

      // Downdate the Gram matrix
      CholUpDown::cholDowndate(L, Eigen::VectorXd::Zero(1), X(all, A_set.head(size_A)), drop_index, l2);

    }
    else
    {

      // If all variables have non-zero coefficients -> return
      if (N == size_A)
      {
        // Return the penalties or the coefficient vector depending on the needs
        if (return_penalties) {
          return pen(seq(0, s)) * (2.0 / l2_sqrt_inv);
        }
        else {
          return beta_hat;
        }
      }
    }
    ++s;
  }

  /* End LARS-EN loop */

  // Return the penalties or the coefficient vector depending on the needs
  if (return_penalties) {
    return pen(seq(0, s)) * (2.0 / l2_sqrt_inv);
  }
  else {
    return beta_hat;
  }
}


#endif /* defined(ENET_SOLVER) */
