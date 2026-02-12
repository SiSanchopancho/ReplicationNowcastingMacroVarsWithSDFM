/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2024-2026 Domenic Franjic
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

#include "Forecast.h"

/* Compute forecast/nowcasts given the higher frequency factor estiamtes */
double Forecast::factorForecaster(const Eigen::MatrixXd& X, const Filtering::KFS_fit& Fit, const int& K, const int& frequency, const int& delay, const int& h, const int& x0_ind)
{
    /* Dummies */

    // Integers
    int tau = 0;

    // Reals
    double x_h = DBL_MAX;

    // Vectors
    Eigen::VectorXd lambda_0 = Eigen::VectorXd::Zero(K), x0_temp = Eigen::VectorXd::Zero((h + 1) / 3);

    // Matrices
    Eigen::MatrixXd F_temp = Eigen::MatrixXd::Zero(K, (h + 1) / 3);;

    /* Forecast according to the frequency of the VOI */
    if (frequency == 4)
    {

        /* Case 1: mixed frequency */ 

        // Aggregate the factor and extract the corresponding quarterly variables if VOI is of qaurterly frequency
        Eigen::MatrixXd F_hat = Eigen::MatrixXd::Zero(K, h + 3);

        // SDFMKFS will loose the first o - 1 observations of the estimated factors, where o is the order of the VAR process in the transition-equation
        // Further, F_hat will need two additional observations at the beginning of the pannel due to Mariano-Murasawa
        F_hat(Eigen::all, Eigen::seq(Fit.order + 1, last)) = Fit.F(Eigen::seq(0, K - 1), Eigen::seq(0, Eigen::last - 1));

        // Aggregate the factors according to Mariano and Murasawa
        // Note: if K == 1 F_hat must be treated as vector for the loop to work.
        // Also: Extract the quarterly obsevrations of the VOI
        // Note: the convention is used that the first observation of the data panel is the first month of the first quarter.
        int t = 4;
        int tt = 0;
        if (K == 1)
        {
            while (t <= h + 2)
            {

                F_temp(tt) = 1. / 3. * (F_hat(t) + 2. * F_hat(t - 1) + 3. * F_hat(t - 2) + 2. * F_hat(t - 3) + F_hat(t - 4));

                x0_temp(tt) = X(t - 2, x0_ind);

                t += 3;
                ++tt;
            }
        }
        else
        {
            while (t <= h + 2)
            {

                F_temp(Eigen::all, tt) = 1. / 3. * (F_hat(Eigen::all, t) + 2. * F_hat(Eigen::all, t - 1) + 3. * F_hat(Eigen::all, t - 2) + 2. * F_hat(Eigen::all, t - 3) + F_hat(Eigen::all, t - 4));

                x0_temp(tt) = X(t - 2, x0_ind);

                t += 3;
                ++tt;
            }
        }

        // Remove the first column since not enough observations are available for aggregation.
        // Further, remove every next observation if the lag order exceeds a multiple of three since this interferes with the estimation of the upcomming warterly observation.
        // Accordingly, decrease the index corresponding to the last observation for which the VOI is observed at forecasting date (tau) and decrease the index of the realisation of the forecast (last_obs).
        tau = static_cast<int>((h + 1) / 3 - std::floor(double(delay) / 3.) - 1);
        for (int o = 0; o < std::floor(double(Fit.order) / 3.) + 1; ++o)
        {
            DataHandle::removeCol(F_temp, 0, false);
            DataHandle::removeElement(x0_temp, 0, false);
            --tau;
        }
        x0_temp.conservativeResize(x0_temp.size() - 1);

        // Save the observation that is aimed to be nowcast
        x_h = X(h, x0_ind);

    }
    else if (frequency == 12)
    {

      throw std::runtime_error("Forecasting for monthly frequencies (12) is currently disabled.");

    }

    // Nowcast the variable of interest
    // Calculate the OLS fit of the loadings of the VOI on the estimated factors
    lambda_0 = ((F_temp(Eigen::all, Eigen::seq(0, tau)) * F_temp(Eigen::all, Eigen::seq(0, tau)).transpose()).llt().solve(Eigen::MatrixXd::Identity(K, K)) * (F_temp(Eigen::all, Eigen::seq(0, tau)) * x0_temp));
    return (x_h - lambda_0.transpose() * F_temp(all, last));
}

/* Compute the rotational distance between the factor estimates and the true factors */
double Forecast::factorDistance(const Eigen::MatrixXd& F_hat, const Eigen::MatrixXd& F)
{
    // Get OLS fit
    Eigen::MatrixXd Beta = (F_hat * F_hat.transpose()).llt().solve(Eigen::MatrixXd::Identity(F_hat.rows(), F_hat.rows())) * (F_hat * F);

    return (F.transpose() - Beta * F_hat).squaredNorm();

}
