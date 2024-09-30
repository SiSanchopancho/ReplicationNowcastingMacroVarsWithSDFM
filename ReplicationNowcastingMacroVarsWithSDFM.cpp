/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2024 Domenic Franjic
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

#define _USE_MATH_DEFINES

 // Defin PI when using g++
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// External includes
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <Eigen>
#include <algorithm>
#include <cctype>
#include <ctime>
#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#endif

// Internal includes
#include "Internals/DataGen.h"
#include "Internals/Developer.h"
#include "Internals/ElNetSolve.h"
#include "Internals/SparsePCA.h"
#include "Internals/DataHandle.h"
#include "Internals/Filtering.h"
#include "Internals/SparseDFM.h"
#include "Internals/CrossVal.h"
#include "Internals/Orders.h"
#include "Internals/Forecast.h"

/* Some helper functions for user interaction */

// Display notices
void displayNotice(const std::string& type) {
    if (type == "c") {
        std::cout << '\n' << "Copying Conditions \n"
            << "This program is free software : you can redistribute itand /or modify it under the "
            << "terms of the GNU General Public License as published by the Free Software Foundation, "
            << "either version 3 of the License, or (at your option) any later version. \n"
            << "This program is distributed in the hope that it will be useful, but WITHOUT "
            << "ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS "
            << "FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details. \n"
            << "You should have received a copy of the GNU General Public License along with "
            << "this program.If not, see https ://www.gnu.org/licenses/." << std::endl << std::endl;
    }
    else if (type == "w") {
        std::cout << '\n' << "15. Disclaimer of Warranty \n"
            << "THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. "
            << "EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND / OR OTHER "
            << "PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED "
            << "OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY "
            << "AND FITNESS FOR A PARTICULAR PURPOSE.THE ENTIRE RISK AS TO THE QUALITY AND "
            << "PERFORMANCE OF THE PROGRAM IS WITH YOU.SHOULD THE PROGRAM PROVE DEFECTIVE, YOU "
            << "ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR, OR CORRECTION." << std::endl << std::endl;
    }
}

// Convert input to lower case letters
std::string toLower(const std::string& input) {
    std::string output = input;
    std::transform(output.begin(), output.end(), output.begin(),
        [](unsigned char c) { return std::tolower(c); });
    return output;
}

// Get user inputs
std::string getUserInput(const std::string& prompt) {
    while (true) {
        std::cout << prompt;
        std::string input;
        std::getline(std::cin, input);
        std::string lowerInput = toLower(input);
        return input;
    }
}

/* Set "hard" parameters */

const int N = 50; // Number of variables
const int T = 3 * (100) + 3; // Number of observations in monthly terms
const int FCH = 3; // Forecasting horizon
const int T_CV = T - FCH; // Cross-Validation set
const int H = static_cast<int>(std::ceil(double(T) / 3.)); // Number of CV nowcasts
const int II = 10; // Number of experiment repetitions
int R_vec[4] = { 1 , 2 , 3 , 4 }; // Array storing the values of factors ought to be checked
double beta_vec[4] = { 4.212120 , 1.737370 , 0.994949 , 0.5 }; // Arrray storing the degrees of correlation ought to be checked
double prob_vec[2] = { 0.0, 0.8 }; // Array storing the degrees of sparsities ought to be checked
const int order = 1; // Order of factor VAR process
bool should_be_corr = true; // Should data be correlated inside the data generating function Yes/No
const bool decorr_errors = true; // Should errors be decorrelated in the SDFM model case?
bool run_once = 0; // Indicator whether the routine has been run once already
const int G_rnd = 500; // Size of the random hyper-parameter search
const int seed = 18092024;
std::mt19937 gen(seed); // RNG

int main() {

    Eigen::initParallel();

    // Detect if stdout is connected to a terminal
    if (isatty(fileno(stdout))) {

        // Output the notice
        std::cout << "ReplicationNowcastingMacroVarsWithSDFM Copyright (C) 2024  Domenic Franjic" << std::endl;
        std::cout << "This program comes with ABSOLUTELY NO WARRANTY; for details type `w'." << std::endl;
        std::cout << "This is free software, and you are welcome to redistribute it" << std::endl;
        std::cout << "under certain conditions; type `c' for details." << std::endl;
        std::string input = toLower(getUserInput("To continue with the program type `continue' or hit enter."));
        while (true) {
            if (input == "w") {
                displayNotice("w");
                input = "";
                input = toLower(getUserInput("To continue with the program type `continue' or hit enter. To see the copying conditions type 'c' "));
            }
            else if (input == "c") {
                displayNotice("c");
                input = "";
                input = toLower(getUserInput("To continue with the program type `continue' or hit enter. To see the disclaimer of waranty type 'w' "));
            }
            else if (input == "" || input == "continue") {
                break;
            }
            else {
                input = toLower(getUserInput("Invalid input. Please enter 'w', 'c' or 'continue'. "));
            }
        }

        // Clear the shell for visibility and re-print the notice
#ifdef _WIN32
        system("cls");
#else
        system("clear");
#endif
        std::cout << "ReplicationNowcastingMacroVarsWithSDFM Copyright (C) 2024  Domenic Franjic" << std::endl;
        std::cout << "This program comes with ABSOLUTELY NO WARRANTY; for details type `w'." << std::endl;
        std::cout << "This is free software, and you are welcome to redistribute it" << std::endl;
        std::cout << "under certain conditions; type `c' for details." << std::endl << std::endl;

    }


    /* Interaction chunk */

    // Initialize values for user input
    int R_last = -1;
    double prob_last = -1, beta_param_last = DBL_MAX;
    int i_start = 0;
    bool rerun = false;

    // Ask wether everything should be re-run
    std::string input = toLower(getUserInput("Do you want to rerun everything? (yes/no): "));
    while (true) {
        if (input == "yes" || input == "y" || input == "true" || input == "1") {
            rerun = true;
            std::cout << std::endl;
            break;
        }
        else if (input == "no" || input == "n" || input == "false" || input == "0") {
            rerun = false;
            std::cout << std::endl;
            break;
        }
        else {
            input = toLower(getUserInput("\nError! Invalid input. Please enter 'yes' or 'no'.: "));
        }
    }


    // If data is not rerun, collect input from the user indicating where the loop should start
    if (!rerun) {
        std::cout << "Where do you want to start?" << std::endl << std::endl;

        // Get the number of factors
        input = getUserInput("Please enter the number of factors: ");
        while (true) {
            try {
                R_last = std::stoi(input);
                if (0 <= R_last) {
                    std::cout << std::endl;
                    break;
                }
                else {
                    input = toLower(getUserInput("\nError! Please enter a non - negative integer: "));
                }
            }
            catch (...) {
                input = toLower(getUserInput("\nError! Invalid input. Please enter an integer value: "));
            }
        }

        // Get the sparsity ratio
        input = getUserInput("Please enter the sparsity ratio (0.0 to 1.0): ");
        while (true) {
            try {
                prob_last = std::stod(input);
                if (0.0 <= prob_last && prob_last <= 1.0) {
                    std::cout << std::endl;
                    break;
                }
                else {
                    input = toLower(getUserInput("\nError! Please enter a value between 0.0 and 1.0: "));
                }
            }
            catch (...) {
                input = toLower(getUserInput("\nInvalid input. Please enter a numerical value.: "));
            }
        }

        // Get the degree of correlation
        input = getUserInput("Please enter the degree of correlation: ");
        while (true) {
            try {
                beta_param_last = std::stod(input);
                if (0.0 <= beta_param_last) {
                    std::cout << std::endl;
                    break;
                }
                else {
                    input = toLower(getUserInput("\nError! Please enter a non-negative: "));
                }
            }
            catch (...) {
                input = toLower(getUserInput("\nInvalid input. Please enter a numerical value.: "));
            }
        }

        // Get the starting position in the loop
        input = getUserInput("Please enter the position in the loop (starting with 0): ");
        while (true) {
            try {
                i_start = std::stoi(input);
                if (i_start >= 0) {
                    std::cout << std::endl;
                    break;
                }
                else {
                    input = toLower(getUserInput("\nError! Please enter a non - negative integer: "));
                }
            }
            catch (...) {
                input = toLower(getUserInput("\nError! Invalid input. Please enter an integer value: "));
            }
        }
    }

    /* Main simulation loop */

    for (int R : R_vec)
    {

        // Skip loops over the number of factors until the approprate starting point is reached
        if (R < R_last)
        {
            continue;
        }

        // Create a subfolder for the results
        std::string res_folder_name = "./results";
        DataHandle::ensureFolderExists(res_folder_name);
        res_folder_name = "./results/results_T_" + std::to_string(T) + "_N_" + std::to_string(N) + "_K_" + std::to_string(R);
        DataHandle::ensureFolderExists(res_folder_name);

        R_last = -1;

        for (double prob : prob_vec)
        {

            // Skip loops over the sparsity ratio until the appropriate starting point is reached
            if (prob < prob_last)
            {
                continue;
            }
            prob_last = -1;

            // Skip loops over the beta parameter until the appropriate starting point is reached
            for (double beta_param : beta_vec)
            {

                if (beta_param_last < beta_param)
                {
                    continue;
                }
                beta_param_last = DBL_MAX;

                /* Initialise dummies for saving the RMSEand forecasting errors of the different models used(ARMA, TP, Dense, Sparse(BIC CV)) */

                // Double matrices
                Eigen::MatrixXd Results = Eigen::MatrixXd::Zero(II, 5);


                // Double vectors
                VectorXd fe_ARMA = VectorXd::Zero(int(std::floor(FCH / 3))), fe_UM = VectorXd::Zero(int(std::floor(FCH / 3))), fe_d = VectorXd::Zero(int(std::floor(FCH / 3))),
                    fe_s_BIC = VectorXd::Zero(int(std::floor(FCH / 3))), fe_s_CV = VectorXd::Zero(int(std::floor(FCH / 3)));

                // Reals
                double MSFE_ARMA = 0.0, MSFE_UM = 0.0, MSFE_d = 0.0, MSFE_s_BIC = 0.0, MSFE_s_CV = 0.0;

                // If any loop has been executed once, set the starting value for the forecast loop to 0
                if (run_once == 1)
                {
                    i_start = 0;
                }

                // When restarting at a later point in the simulation loop load the old results
                if (i_start != 0)
                {

                    std::string in_name = DataHandle::namer(res_folder_name + "/RMSFE", R, N, T, prob, beta_param, should_be_corr, FCH, T_CV, G_rnd, seed);
                    Results = DataHandle::dataLoadCSV(in_name + ".txt");

                }

                /* Main repition loop for each parametrisation */

                for (int i = i_start; i < II; ++i)
                {

                    // Print the Current run number and parametrisation for tracking
                    std::cout << '\n' << "############## Run number " << i << " for beta = " << beta_param << ", prob = " << prob << ", and R = " << R << " ############" << '\n';

                    // Object for storing the simulated Factor Models
                    DataGen::FM SFM(R, T, N);

                    /* Simulating the data */

                    {
                        /* Create parameters inside this scope for memory management */
                        int burn_in = 250; // Rank of the state error
                        bool quarterfy = true; // Quarterfy some of the data Yes/No
                        double m = 1. / double(N); // Proportion of data that is ought to be quarterfied
                        Eigen::MatrixXd S = Eigen::MatrixXd::Identity(R, R); // Varaince-Covariance matrix of the factor VAR process errors
                        VectorXd mu_e = VectorXd::Zero(N); // Factor means
                        Eigen::MatrixXd Lambda_sim = DataGen::MVData(gen, N, R); // Factor loadings matrix
                        Eigen::MatrixXd Sigma_e_sim = Eigen::MatrixXd::Identity(N, N); // Measurement error variance-covariance matrix
                        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(R, R * order); // Transition process VAR matrix
                        for (int o = 0; o < order; ++o)
                        {
                            A.block(0, R * o, R, R) = (1. / (2. + double(o))) * Eigen::MatrixXd::Identity(R, R);
                            A.block(0, R * o, R, R).diagonal() *= pow(-1., o);
                        }

                        // If the population loadings matrix is ought to be sparse, set zeros and rescale the columns (currently not in use but kept for potential later usage)
                        if (R == 2)
                        {
                            if (prob != 0.)
                            {

                                Lambda_sim.topRightCorner(std::floor((double(N) * prob)), 1).setZero();
                                Lambda_sim.bottomLeftCorner(std::floor((double(N) * prob)), 1).setZero();

                                for (int k = 0; k < R; ++k)
                                {
                                    Lambda_sim.col(k).normalize();
                                    Lambda_sim.col(k) *= std::sqrt(double(N));
                                }

                                Lambda_sim.topRightCorner(std::floor((double(N) * prob)), 1).setZero();
                                Lambda_sim.bottomLeftCorner(std::floor((double(N) * prob)), 1).setZero();

                            }
                        }
                        else if (R == 1)
                        {
                            if (prob != 0.)
                            {

                                Lambda_sim(seq(0, std::floor(double(N) * prob)), 0).setZero();

                                Lambda_sim.col(0).normalize();
                                Lambda_sim.col(0) *= std::sqrt(double(N));

                                Lambda_sim(seq(0, std::floor(double(N) * prob)), 0).setZero();

                            }
                        }

                        // Simulate the DFM
                        DataGen::staticFM(SFM, T, N, S, Lambda_sim, mu_e, Sigma_e_sim, A, gen, order, quarterfy, should_be_corr, beta_param, m, R, burn_in, true);

                    }

                    // Standardise the data
                    DataHandle::standardise(SFM.X);

                    /* Forecast preperations */

                    // Create a frequency indicator
                    int  months = int((SFM.frequency.array() == 12).count()), ind_fre_indicator = 0, max_order = 10;
                    Eigen::VectorXi  freq_ind = Eigen::VectorXi::Zero(months);
                    for (int n = 0; n < N; ++n)
                    {
                        if (SFM.frequency(n) == 12)
                        {
                            freq_ind(ind_fre_indicator) = n;
                            ++ind_fre_indicator;
                        }
                    }

                    // Create an artificial delay vector (in months) where VOI is lagged by 3 units (1 quarter)
                    Eigen::VectorXi delay = Eigen::VectorXi::Zero(N);
                    delay(0) = 3;


                    /* Validation */

                    // Create results object
                    CrossVal::CV_res CVs_BIC(R), CVs_CV(R);

                    // Cross-validation
                    CVs_CV = CrossVal::parallelRandCrossVal(0, SFM.X(seq(0, T_CV), all), H, R, G_rnd, delay, SFM.date(seq(0, T_CV)), SFM.frequency, gen, decorr_errors, 10,
                        -6, +6, 5, 10000000, "BIC", 10e-15, 10e-6, true, true);

                    // Information criterion
                    CVs_BIC = CrossVal::parallelRandBIC(SFM.X(seq(0, T_CV), all), R, G_rnd, delay, SFM.date(seq(0, T_CV)), SFM.frequency, gen, decorr_errors, 10, -6, +6, 5, 10000000, "BIC", 10e-15,
                        10e-6, true, true);

                    /* Forecast */

                    /* Forecasting Dummies */

                    // Integers
                    int VOI_ind = 0, number_of_nforecast_indicator = 0, max_iterations = 1000, running = 0;

                    // Matrices
                    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(R, R);

                    // Double vectors
                    VectorXd lambda_0_s_BIC = VectorXd::Zero(R), lambda_0_s_CV = VectorXd::Zero(R), lambda_0_d = VectorXd::Zero(R), lambda_0_TP = VectorXd::Zero(R),
                        l1 = VectorXd::Constant(R, NAN);

                    // Integer vectors
                    Eigen::VectorXi curr_sel = Eigen::VectorXi::Zero(R), ind_vec = Eigen::VectorXi::Zero(R), fc_dates = Eigen::VectorXi::Zero(H);

                    // Create a vector of "forecasting dates" depending on the frequency of the VOI
                    if (SFM.frequency(VOI_ind) == 4)
                    {

                        fc_dates.resize(int(std::round(double(FCH) / 3.)));
                        int tt = 0;
                        for (int t : SFM.date(seq(T - FCH + 1, T - 1)))
                        {
                            if ((t + 1) % 3 == 0)
                            {
                                fc_dates(tt) = t;
                                ++tt;
                            }
                        }
                    }
                    else if (SFM.frequency(VOI_ind) == 12)
                    {
                        fc_dates = Eigen::VectorXi::LinSpaced(FCH, T - FCH, T - 1);
                    }

                    /* Forecasting */

                    std::cout << '\n' << "Forecasting. Please stand by." << "      ";

                    /* Main forecast loop */

                    for (int h : fc_dates)
                    {

                        // Indicator showing that the program is doing something
                        if (running % 4 == 0)
                        {
                            std::cout << "\b\b" << " -" << std::flush;
                        }
                        else if (running % 4 == 1)
                        {
                            std::cout << "\b\b" << " \\" << std::flush;
                        }
                        else if (running % 4 == 2)
                        {
                            std::cout << "\b\b" << " |" << std::flush;
                        }
                        else if (running % 4 == 3)
                        {
                            std::cout << "\b\b" << " /" << std::flush;
                        }
                        ++running;

                        // Estimate the dense model using only the monthly observations up to T = h
                        Filtering::KFS_fit DFM_d;
                        SparseDFM::DFMKFS(DFM_d, SFM.X(seq(0, h), freq_ind), delay, R, 10);

                        // Estimate the sparse model validated via BIC
                        Filtering::KFS_fit DFM_s_BIC;
                        SparseDFM::SDFMKFS(DFM_s_BIC, SFM.X(seq(0, h), freq_ind), delay(freq_ind), CVs_BIC.sel, R, max_order, decorr_errors, "BIC", "LARS", CVs_BIC.l2, l1, 0., NAN, 0.01, INT_MIN, max_iterations, INT_MIN, 10e-15,
                            0, 0.0001, 10e-7, 0);

                        // Estimate the sparse model validated via CV
                        Filtering::KFS_fit DFM_s_CV;
                        SparseDFM::SDFMKFS(DFM_s_CV, SFM.X(seq(0, h), freq_ind), delay(freq_ind), CVs_CV.sel, R, max_order, decorr_errors, "BIC", "LARS", CVs_CV.l2, l1, 0., NAN, 0.01, INT_MIN, max_iterations, INT_MIN, 10e-15,
                            0, 0.0001, 10e-7, 0);

                        // Nowcacsting
                        fe_s_BIC(number_of_nforecast_indicator) = Forecast::factorForecaster(SFM.X, DFM_s_BIC, R, SFM.frequency(VOI_ind), delay(VOI_ind), h, VOI_ind);
                        fe_s_CV(number_of_nforecast_indicator) = Forecast::factorForecaster(SFM.X, DFM_s_CV, R, SFM.frequency(VOI_ind), delay(VOI_ind), h, VOI_ind);
                        fe_d(number_of_nforecast_indicator) = Forecast::factorForecaster(SFM.X, DFM_d, R, SFM.frequency(VOI_ind), delay(VOI_ind), h, VOI_ind);

                        /* Estimate some benchmarks */

                        // Forecast using an AR(P) model
                        Eigen::MatrixXd x_ARMA = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(std::floor(double(h - delay.maxCoeff() + 1) / 3.)), 1);
                        int last_obs_ARMA = static_cast<int>(x_ARMA.rows());
                        int tt = 0;
                        for (int ttt = 0; ttt <= h - delay.maxCoeff(); ++ttt)
                        {

                            if ((ttt + 1) % 3 == 0)
                            {
                                x_ARMA(tt, 0) = SFM.X(ttt, 0);
                                ++tt;
                            }

                        }

                        // Estimate the order
                        int order_var = VARorder<BIC>(x_ARMA.transpose(), 10, 10e-15);

                        // Calculate the VAR coefficient matrix Phi
                        Eigen::MatrixXd X_l(last_obs_ARMA, order_var), Phi = Eigen::MatrixXd::Zero(order_var, order_var);

                        for (int oo = 0; oo < order_var; ++oo)
                        {

                            X_l.block(oo, oo, last_obs_ARMA - oo, 1) = x_ARMA.block(0, 0, last_obs_ARMA - oo, 1);

                        }
                        for (int oo = 1; oo < order_var; ++oo)
                        {

                            DataHandle::removeRow(X_l, 0, 0);
                            DataHandle::removeRow(x_ARMA, 0);

                        }

                        // Calculate estimator of the factor VARMA process coefficient matrix Phi
                        Eigen::MatrixXd X_t = X_l(seq(1, last), all);
                        Eigen::MatrixXd X_t_lag = X_l(seq(0, last - 1), all);
                        Eigen::MatrixXd Idento = Eigen::MatrixXd::Identity(order_var, order_var);
                        Phi = ((X_t_lag.transpose() * X_t_lag).llt().solve(Idento) * X_t_lag.transpose() * X_t).transpose();
                        Phi.bottomLeftCorner((order_var - 1), (order_var - 1)) = Eigen::MatrixXd::Identity((order_var - 1), (order_var - 1));
                        Phi.bottomRightCorner((order_var - 1), 1) = Eigen::MatrixXd::Zero((order_var - 1), 1);

                        // Forecast
                        Eigen::MatrixXd x_hat_AR = X_l.bottomLeftCorner(1, order_var) * Phi.row(0).transpose();
                        fe_ARMA(number_of_nforecast_indicator) = SFM.X(h, VOI_ind) - x_hat_AR(0, 0);


                        // Unconditional mean forecast
                        fe_UM(number_of_nforecast_indicator) = SFM.X(h, VOI_ind) - x_ARMA.mean();

                        ++number_of_nforecast_indicator;

                    }

                    /* End main forecast loop */

                    // Print and save the results

                    MSFE_d = ((fe_d.array() * fe_d.array()).mean());
                    MSFE_s_BIC = ((fe_s_BIC.array() * fe_s_BIC.array()).mean());
                    MSFE_s_CV = ((fe_s_CV.array() * fe_s_CV.array()).mean());
                    MSFE_ARMA = ((fe_ARMA.array() * fe_ARMA.array()).mean());
                    MSFE_UM = ((fe_UM.array() * fe_UM.array()).mean());

                    // Store the final RMSFE results for each run
                    Results(i, 0) = MSFE_d;
                    Results(i, 1) = MSFE_s_BIC;
                    Results(i, 2) = MSFE_s_CV;
                    Results(i, 3) = MSFE_ARMA;
                    Results(i, 4) = MSFE_UM;

                    // Safe the results externally
#ifdef _MSC_VER // Check if the compiler is MSVC
#pragma warning(disable : 4996) // Disable MSVC warning due to strcpy
#endif
                    std::string out_name = DataHandle::namer(res_folder_name + "/RMSFE", R, N, T, prob, beta_param, should_be_corr, FCH, T_CV, G_rnd, seed);
                    DataHandle::dataSaveCSV(Results, out_name, ".csv");

                    // Reset the storing variables to zero for the next run just to be sure
                    fe_ARMA.setZero();  fe_UM.setZero(); fe_d.setZero(); fe_s_BIC.setZero(); fe_s_CV.setZero();
                    MSFE_d = 0.0; MSFE_s_BIC = 0.0; MSFE_s_CV = 0.0; MSFE_ARMA = 0.0; MSFE_UM = 0.0;

                }

                /* End main repition loop for each parametrisation */

                // Print the results for this parametrisation
                std::cout << '\n';
                std::cout << '\n' << "The MSNE for the Dense model is: " << '\t' << Results.col(0).mean() << '\n'
                    << "The MSNE for the Sparse BIC model is: " << '\t' << Results.col(1).mean() << '\n'
                    << "The MSNE for the Sparse CV model is: " << '\t' << Results.col(2).mean() << '\n'
                    << "The MSNE for the ARMA model is: " << '\t' << Results.col(3).mean() << '\n'
                    << "The MSNE for the UM model is: " << '\t' << Results.col(4).mean() << '\n';

                // Indicate that the loop has run at least once
                run_once = 1;

            }
        }
    }

    /* End main simulation loop */

    // After the simulation has run, make sure that the shell is not closing immediately
    std::cout << "Enough? If so press any key.";
    std::cin.ignore();
    std::cin.get();

    return 0;
}

