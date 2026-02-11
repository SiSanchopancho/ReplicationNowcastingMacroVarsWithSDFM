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

#include "CrossVal.h"

/* Random Time Series Cross-Validation */
CrossVal::CV_res CrossVal::randCrossVal(
    const int& x0_ind, // Index of the variable of interest
    const Eigen::MatrixXd& X, // (T x N) Data matrix
    const int& H_in, // Number of forecasts for the cross-validation
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum number of models ought to be checked
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const int& sel_min, // Lower bound of the distribution of the number of selected variables
    int sel_max, // Upper bound of the distribution of the number of selected variables
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Do you want to time cross-validation (Y/N)
    const bool& log // Talk to me
)
{

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    /* Initialisation */

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count()), H = H_in;

    if (months < sel_max) {

        sel_max = months;

    }

    if (H % 3 != 0) {

        while (H % 3 != 0) {

            H += 1;

        }

    }

    // Matrices
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(R, R);

    // Vectors
    Eigen::VectorXi freq_ind = Eigen::VectorXi::Zero(months), fc_dates = Eigen::VectorXi::Zero(H);

    // Fits
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    // Check for miss-handling

    // Check whether the first observation of the panel is the first month of the quarter
    if (frequency(x0_ind) == 4 && (X(0, x0_ind) != X(1, x0_ind) || X(0, x0_ind) != X(2, x0_ind)))
    {

        std::cout << '\n' << "Error! One of two things went wrong!" << '\n'
            << "If the variable of interest is of quarterly frequency, the first observation of the panel must be the first month of the first quarter. Try removing some observatinos and restart the calculations." << '\n'
            << "Otherwise, the convention is used that the missing observations for the first and second month of the qaurterly variables are imputed by the corresponding value of the quarter observed at the last month of the quarter. Check whether this is the case." << '\n';

        return results;

    }

    // Create a vector indicating which variables are of the lowest frequency (aon only monthly)

    for (int n = 0; n < N; ++n)
    {
        if (frequency(n) == 12)
        {
            freq_ind(ind) = n;
            ++ind;
        }
    }

    // Create pseudo date vectors based on the frequency of VOI

    if (frequency(x0_ind) == 4)
    {

        fc_dates.resize(static_cast<Eigen::Index>(std::floor(double(H) / 3.)));

        int tt = 0;

        for (int t : date(seq(T - H + 1, T - 1)))
        {
            if ((t + 1) % 3 == 0)
            {

                fc_dates(tt) = t;
                ++tt;

            }
        }
    }
    else if (frequency(x0_ind) == 12)
    {

        fc_dates = Eigen::VectorXi::LinSpaced(H, T - H, T - 1);

    }

    // Print the fc dates for monitoring purposes

    int cv_dates = static_cast<int>(fc_dates.size());
    Eigen::MatrixXi Check_dates = Eigen::MatrixXi::Zero(int(floor(double(cv_dates) / 2.)), 2);
    Check_dates.col(0) = fc_dates.head(int(floor(double(cv_dates) / 2.)));
    Check_dates.col(1) = fc_dates.tail(int(floor(double(cv_dates) / 2.)));

    // Create dummies to store the fes and squares of the residuals

    Eigen::MatrixXd fes = Eigen::MatrixXd::Zero(fc_dates.size(), 2), all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::VectorXd res_sq = Eigen::VectorXd::Zero(fc_dates.size());
    Eigen::MatrixXi all_sel = Eigen::MatrixXi::Constant(R, max_models, months);
    Eigen::VectorXi curr_sel_local = Eigen::VectorXi::Zero(R);

    /* Main cv loop */ 

    if (log) { std::wcout << "\n" << "Currently cross-validating. Please stand by." << "\n" << "\rProgress: " << 0 << "% " << std::flush; }

    // Initialize the generators for the random hyper-parameters
    std::uniform_real_distribution<double> dis_l2(l2_min, l2_max);
    std::uniform_int_distribution<int> dis_l1(sel_min, sel_max);


    for (int i = 1; i < max_models; ++i) {
        all_l2(i) = std::pow(10, dis_l2(gen));
        for (int k = 0; k < R; ++k) {
            all_sel(k, i) = dis_l1(gen);
        }
    }

    for (int i = 0; i < max_models; ++i) {

        double l2 = all_l2(i);
        Eigen::VectorXi curr_sel_local = all_sel.col(i);

        /* Calculate the state - space model for the current hyperparameters */

        int ind2 = 0;
        for (int h : fc_dates)
        {

            // Estimate the SDFM
            Filtering::KFS_fit DFM;
            SparseDFM::SDFMKFS(DFM, X(seq(0, h), freq_ind), delay(freq_ind), curr_sel_local, R, order, decorr_errors, crit, "LARS", l2, Eigen::VectorXd::Constant(R, NAN), 0., NAN, NAN, 1, 1000, INT_MIN, comp_null,
                0, 0.0001, conv_threshold, 0);

            // Store the estimates according to the frequency of the VOI
            fes(ind2, 1) = Forecast::factorForecaster(X, DFM, R, frequency(x0_ind), delay(x0_ind), h, x0_ind);
            ++ind2;
            std::cout << "\b\b" << std::flush;
        }

        // Save the residuals squares
        res_sq = (fes.col(1).array() * fes.col(1).array()).matrix();

        // Save the results of this run if it is the mean of the squared residuals is smaller then the last best run. Else, discard the results.
        double new_min_mse = res_sq.mean();

        if (new_min_mse < results.min_mse)
        {
            results.l2 = l2;
            results.sel = curr_sel_local;
            results.min_mse = new_min_mse;
            fes.col(0) = fes.col(1);
        }

        // Increment the progress counter
        if (log) { std::cout << "\rProgress: " << (100 * i / max_models) << "% " << std::flush; }

    }

    /* End of cv loop */

    if (log) { std::cout << "\rDone. The following parametrisation is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{CV Error} = " << results.min_mse << '\n'; }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    if (timer)
    {

        // Print the duration
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;

    }

    return results;
}

/* Parallelised Random Time Series Cross-Validation */
CrossVal::CV_res CrossVal::parallelRandCrossVal(
    const int& x0_ind, // Index of the variable of interest
    const Eigen::MatrixXd& X, // (T x N) Data matrix
    const int& H_in, // Number of forecasts for the cross-validation
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum number of models ought to be checked
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const int& sel_min, // Lower bound of the distribution of the number of selected variables
    int sel_max, // Upper bound of the distribution of the number of selected variables
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Do you want to time cross-validation (Y/N)
    const bool& log // Talk to me
)
{

    // Parallel rnd CV
    auto start = std::chrono::high_resolution_clock::now();

    /* Dummies */

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count()), H = H_in;

    if (months < sel_max) {
        sel_max = months;
    }

    if (H % 3 != 0) {

        while (H % 3 != 0) {

            H += 1;

        }

    }

    // Matrices and vectors
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(R, R);
    Eigen::VectorXi freq_ind = Eigen::VectorXi::Zero(months), fc_dates = Eigen::VectorXi::Zero(H);

    // Result obj.
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    // Funciton misshandling
    if (frequency(x0_ind) == 4 && (X(0, x0_ind) != X(1, x0_ind) || X(0, x0_ind) != X(2, x0_ind))) {
        std::cout << '\n' << "Error! One of two things went wrong!" << '\n'
            << "If the variable of interest is of quarterly frequency, the first observation of the panel must be the first month of the first quarter. Try removing some observatinos and restart the calculations." << '\n'
            << "Otherwise, the convention is used that the missing observations for the first and second month of the qaurterly variables are imputed by the corresponding value of the quarter observed at the last month of the quarter. Check whether this is the case." << '\n';

        return results;
    }

    // Create pseudo frequency indicator vector
    for (int n = 0; n < N; ++n) {
        if (frequency(n) == 12) {
            freq_ind(ind) = n;
            ++ind;
        }
    }
    if (frequency(x0_ind) == 4) {
        fc_dates.resize(static_cast<Eigen::Index>(std::floor(double(H) / 3.)));
        int tt = 0;
        for (int t : date(seq(T - H + 1, T - 1))) {
            if ((t + 1) % 3 == 0) {
                fc_dates(tt) = t;
                ++tt;
            }
        }
    }
    else if (frequency(x0_ind) == 12) {
        fc_dates = Eigen::VectorXi::LinSpaced(H, T - H, T - 1);
    }

    // Set up the set of random models
    std::uniform_real_distribution<> dis_l2(l2_min, l2_max);
    std::uniform_int_distribution<> dis_l1(sel_min, sel_max);

    Eigen::VectorXd all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::MatrixXi all_sel = Eigen::MatrixXi::Constant(R, max_models, months);


    for (int i = 1; i < max_models; ++i) {
        all_l2(i) = std::pow(10, dis_l2(gen));
        for (int k = 0; k < R; ++k) {
            all_sel(k, i) = dis_l1(gen);
        }
    }

    /* cv loop */

    if (log) { std::cout << "\nCurrently cross-validating. Please stand by." << std::endl; }

    // Serialise Eigen operations within the parallelisation according to https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
    Eigen::setNbThreads(1);

    // Parallel looping
#pragma omp parallel for schedule(dynamic) shared(results)
    for (int i = 0; i < max_models; ++i) {

        // Initialise these objects for each thread locally
        double l2 = all_l2(i);
        Eigen::VectorXi curr_sel_local = all_sel.col(i);
        Eigen::MatrixXd fes_local = Eigen::MatrixXd::Zero(fc_dates.size(), 2);
        Eigen::VectorXd res_sq_local = Eigen::VectorXd::Zero(fc_dates.size());

        // Compute the SDFM
        int ind2 = 0;
        for (int h : fc_dates) {
            Filtering::KFS_fit DFM;
            SparseDFM::SDFMKFS(DFM, X(seq(0, h), freq_ind), delay(freq_ind), curr_sel_local, R, order, decorr_errors, crit, "LARS", l2, Eigen::VectorXd::Constant(R, NAN), 0., NAN, NAN, 1, 1000, INT_MIN, comp_null, 0, 0.0001, conv_threshold, 0);
            fes_local(ind2, 1) = Forecast::factorForecaster(X, DFM, R, frequency(x0_ind), delay(x0_ind), h, x0_ind);
            ++ind2;
        }

        // Store the results
        res_sq_local = (fes_local.col(1).array() * fes_local.col(1).array()).matrix();

        double mse_local = res_sq_local.mean();

        // Compare results between threads in critical environment
#pragma omp critical
        {
            if (mse_local < results.min_mse) {
                results.l2 = l2;
                results.sel = curr_sel_local;
                results.min_mse = mse_local;
                fes_local.col(0) = fes_local.col(1);
            }
        }
    }

    /* end cv loop */

    // Enable Eigen's internal parallelisation again
    Eigen::setNbThreads(0);

    if (log) { std::cout << "\rDone. The following parametrisation is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{CV Error} = " << results.min_mse << "\n"; }

    // Print duration if needed
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    if (timer) {
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;
    }

    return results;
}



CrossVal::CV_res CrossVal::randBIC(
    const Eigen::MatrixXd& X_in, // (T x N) Data matrix
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum Number of models ought to be tested
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const int& sel_min, // Lower bound of the distribution of the number of selected variables
    int sel_max, // Upper bound of the distribution of the number of selected variables
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Timing the duration of the validaiton loop
    const bool& log // Talk to me 
)
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    /* Dummies */ 

    // Matrices
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(R), static_cast<Eigen::Index>(R)), X = X_in;

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count()), progress = 0;

    // Vectors
    Eigen::VectorXd lambda_0 = Eigen::VectorXd::Zero(R);
    Eigen::VectorXi ind_vec = Eigen::VectorXi::Zero(R), freq_ind = Eigen::VectorXi::Zero(months), curr_sel = Eigen::VectorXi::Zero(R);

    // Fits
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    if (months < sel_max) {

        sel_max = months;

    }

    // Create a vector indicating which variables are of the lowest frequency (aon only monthly)
    for (int n = 0; n < N; ++n)
    {
        if (frequency(n) == 12)
        {
            freq_ind(ind) = n;
            ++ind;
        }
    }

    // Initialize the random model sets
    std::uniform_real_distribution<> dis_l2(l2_min, l2_max);
    std::uniform_int_distribution<> dis_l1(sel_min, sel_max);


    Eigen::VectorXd all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::MatrixXi all_sel = Eigen::MatrixXi::Constant(R, max_models, months);

    for (int i = 1; i < max_models; ++i) {
        all_l2(i) = std::pow(10, dis_l2(gen));
        for (int k = 0; k < R; ++k) {
            all_sel(k, i) = dis_l1(gen);
        }
    }

    /* Validation Loop */

    for (int i = 0; i < max_models; ++i) {

        double l2 = all_l2(i);
        curr_sel = all_sel.col(i);

        // Calculate the state-space model for the current hyperparameters
        // Internal constants
        const int max_iterations = 1000;

        // Reset the predictor matrix in case the decorrelating transformation has been used in the last step
        if (decorr_errors)
        {
            X = X_in;
        }

        // Compute the SDFM
        Filtering::KFS_fit DFM;
        SparseDFM::SDFMKFS(DFM, X(all, freq_ind), delay(freq_ind), curr_sel, R, order, decorr_errors, crit, "LARS", l2, Eigen::VectorXd::Constant(R, NAN), 0., NAN, NAN, 1, max_iterations, INT_MIN, comp_null,
            0, 0.0001, conv_threshold, 0);

        if (decorr_errors)
        {
            X(all, freq_ind) = (X(all, freq_ind) * DFM.C.transpose()).eval();
        }

        /* Compute the BIC */
        double BIC_curr = std::log((1. / (double(months) * double(T)) * (X(Eigen::seq(DFM.order - 1, last), freq_ind) - (DFM.Lambda_hat * DFM.F(all, Eigen::seq(0, last - 1))).transpose()).squaredNorm())) + double((DFM.Lambda_hat.array() != 0).count()) * (std::log(double(months) * double(T))) / (double(months) * double(T));

        // Save the results of this run if it is the mean of the squared residuals is smaller than the last best run. Else, discard the results.
        if (BIC_curr < results.min_mse)
        {
            results.l2 = l2;
            results.sel = curr_sel;
            results.min_mse = BIC_curr;
        }


        // Increment the shared progress counter
        progress++;

        if (log) { std::cout << "\rProgress: " << (100 * i / max_models) << "% " << std::flush; }
    }

    /* End Validation Loop */

    if (log) { std::cout << "\rDone. The following parametrization is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{BIC} = " << results.min_mse << '\n'; }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Print the duration if wished for
    if (timer)
    {
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;
    }

    return results;
}

CrossVal::CV_res CrossVal::parallelRandBIC(
    const Eigen::MatrixXd& X_in, // (T x N) Data matrix
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum number of models ought to be tested
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const int& sel_min, // Lower bound of the distribution of the number of selected variables
    int sel_max, // Upper bound of the distribution of the number of selected variables
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Time the validation loop
    const bool& log // Talk to me
)
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    /* Dummies */

    // Matrices
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(R), static_cast<Eigen::Index>(R)), X = X_in;

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count());

    // Vectors
    Eigen::VectorXd lambda_0 = Eigen::VectorXd::Zero(R);
    Eigen::VectorXi freq_ind = Eigen::VectorXi::Zero(months);

    // Fits
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    if (months < sel_max) {

        sel_max = months;

    }

    // Create a vector indicating which variables are of the lowest frequency (aon only monthly)
    for (int n = 0; n < N; ++n)
    {
        if (frequency(n) == 12)
        {
            freq_ind(ind) = n;
            ++ind;
        }
    }

    // Initialize the random model sets
    std::uniform_real_distribution<> dis_l2(l2_min, l2_max);
    std::uniform_int_distribution<> dis_l1(sel_min, sel_max);

    Eigen::VectorXd all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::MatrixXi all_sel = Eigen::MatrixXi::Constant(R, max_models, months);

    for (int i = 1; i < max_models; ++i) {
        all_l2(i) = std::pow(10, dis_l2(gen));
        for (int k = 0; k < R; ++k) {
            all_sel(k, i) = dis_l1(gen);
        }
    }

    /* Validation loop */

    if (log) { std::cout << "\nCurrently cross-validating. Please stand by." << std::endl; }

    // Serialise Eigen operations within the parallelisation according to https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
    Eigen::setNbThreads(1);

    // Parallel looping
#pragma omp parallel for schedule(dynamic) shared(results)
    for (int i = 0; i < max_models; ++i) {

        // Set local thread variables
        const double l2 = all_l2(i);
        const Eigen::VectorXi curr_sel_local = all_sel.col(i);

        // Reset the predictor matrix in case the decorrelating transformation has been used in the last step (locally!)
        Eigen::MatrixXd X_local = X_in;

        // Compute the SDFM
        Filtering::KFS_fit DFM;
        SparseDFM::SDFMKFS(DFM, X_local(all, freq_ind), delay(freq_ind), curr_sel_local, R, order, decorr_errors, crit, "LARS", l2, Eigen::VectorXd::Constant(R, NAN), 0., NAN, NAN, 1, 1000, INT_MIN, comp_null,
            0, 0.0001, conv_threshold, 0);

        if (decorr_errors)
        {
            X_local(all, freq_ind) = (X_local(all, freq_ind) * DFM.C.transpose()).eval();

        }

        // Compute the BIC
        double BIC_curr = std::log((1. / (double(months) * double(T)) * (X_local(Eigen::seq(DFM.order - 1, last), freq_ind) - (DFM.Lambda_hat * DFM.F(all, Eigen::seq(0, last - 1))).transpose()).squaredNorm())) + double((DFM.Lambda_hat.array() != 0).count()) * (std::log(double(months) * double(T))) / (double(months) * double(T));

        // Compare results between threads in critical environment
#pragma omp critical
        {
            if (BIC_curr < results.min_mse)
            {
                results.l2 = l2;
                results.sel = curr_sel_local;
                results.min_mse = BIC_curr;
            }
        }
    }

    /* End validaiton loop */

    // Enable Eigen's internal parallelisation again
    Eigen::setNbThreads(0);

    // End of the loop
    if (log) { std::cout << "\rDone. The following parametrization is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{BIC} = " << results.min_mse << '\n'; }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Print the duration if wished
    if (timer)
    {
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;
    }

    return results;
}

CrossVal::CV_res CrossVal::adjPenaltyParallelRandBIC(
    const Eigen::MatrixXd& X_in, // (T x N) Data matrix
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum number of models ought to be tested
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const int& sel_min, // Lower bound of the distribution of the number of selected variables
    int sel_max, // Upper bound of the distribution of the number of selected variables
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Time the validation loop
    const bool& log // Talk to me
)
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    /* Dummies */

    // Matrices
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(R), static_cast<Eigen::Index>(R)), X = X_in;

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count());

    // Vectors
    Eigen::VectorXd lambda_0 = Eigen::VectorXd::Zero(R);
    Eigen::VectorXi freq_ind = Eigen::VectorXi::Zero(months);

    // Fits
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    if (months < sel_max) {

        sel_max = months;

    }

    // Create a vector indicating which variables are of the lowest frequency (aon only monthly)
    for (int n = 0; n < N; ++n)
    {
        if (frequency(n) == 12)
        {
            freq_ind(ind) = n;
            ++ind;
        }
    }

    // Initialize the random model sets
    std::uniform_real_distribution<> dis_l2(l2_min, l2_max);
    std::uniform_int_distribution<> dis_l1(sel_min, sel_max);

    Eigen::VectorXd all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::MatrixXi all_sel = Eigen::MatrixXi::Constant(R, max_models, months);

    for (int i = 1; i < max_models; ++i) {
        all_l2(i) = std::pow(10, dis_l2(gen));
        for (int k = 0; k < R; ++k) {
            all_sel(k, i) = dis_l1(gen);
        }
    }

    /* Validation loop */

    if (log) { std::cout << "\nCurrently cross-validating. Please stand by." << std::endl; }

    // Serialise Eigen operations within the parallelisation according to https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
    Eigen::setNbThreads(1);

    // Parallel looping
#pragma omp parallel for schedule(dynamic) shared(results)
    for (int i = 0; i < max_models; ++i) {

        // Set local thread variables
        const double l2 = all_l2(i);
        const Eigen::VectorXi curr_sel_local = all_sel.col(i);

        // Reset the predictor matrix in case the decorrelating transformation has been used in the last step (locally!)
        Eigen::MatrixXd X_local = X_in;

        // Compute the SDFM
        Filtering::KFS_fit DFM;
        SparseDFM::SDFMKFS(DFM, X_local(all, freq_ind), delay(freq_ind), curr_sel_local, R, order, decorr_errors, crit, "LARS", l2, Eigen::VectorXd::Constant(R, NAN), 0., NAN, NAN, 1, 1000, INT_MIN, comp_null,
            0, 0.0001, conv_threshold, 0);

        if (decorr_errors)
        {
            X_local(all, freq_ind) = (X_local(all, freq_ind) * DFM.C.transpose()).eval();

        }

        // Compute the 
        double sqrt_NT = std::sqrt(double(months) * double(T));
        double BIC_curr = std::log((1. / (double(months) * double(T)) * (X_local(Eigen::seq(DFM.order - 1, last), freq_ind) - (DFM.Lambda_hat * DFM.F(all, Eigen::seq(0, last - 1))).transpose()).squaredNorm())) + double((DFM.Lambda_hat.array() != 0).count()) * (std::log(sqrt_NT)) / (sqrt_NT);

        // Compare results between threads in critical environment
#pragma omp critical
        {
            if (BIC_curr < results.min_mse)
            {
                results.l2 = l2;
                results.sel = curr_sel_local;
                results.min_mse = BIC_curr;
            }
        }
    }

    /* End validaiton loop */

    // Enable Eigen's internal parallelisation again
    Eigen::setNbThreads(0);

    // End of the loop
    if (log) { std::cout << "\rDone. The following parametrization is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{BIC} = " << results.min_mse << '\n'; }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Print the duration if wished
    if (timer)
    {
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;
    }

    return results;
}

CrossVal::CV_res CrossVal::parallelRandSinglePenaltyBIC(
    const Eigen::MatrixXd& X_in, // (T x N) Data matrix
    const int& R, // Number of factors of the DFm
    const int& max_models, // Maximum number of models ought to be tested
    const Eigen::VectorXi& delay, // Vector of delays
    const Eigen::VectorXi& date, // date vector
    const Eigen::VectorXi& frequency, // frequency vector
    std::mt19937& gen, // RNG
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const double& l1_max, // Maximum possible l1 penalty
    const double& l1_min, // Minimum possible l1 penalty
    const double& l2_min, // Lower bound of the l2 penalty exponent distribution
    const double& l2_max, // Upper bound of the l2 penalty exponent distribution
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer, // Timein the main validation loop
    const bool& log // Talk to me
)
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    /* Dummies */
    
    // Matrices
    Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(R), static_cast<Eigen::Index>(R)), X = X_in, all_l1 = Eigen::MatrixXd::Constant(static_cast<Eigen::Index>(R), static_cast<Eigen::Index>(max_models), 0.0);

    // Integers
    int ind = 0, T = static_cast<int>(X.rows()), N = static_cast<int>(X.cols()), months = static_cast<int>((frequency.array() == 12).count());

    // Reals
    double rnd_li_given_l2 = 0.0;

    // Booleans
    bool l1_max_nan = std::isnan(l1_max), l1_min_nan = std::isnan(l1_min);

    // Vectors
    Eigen::VectorXd lambda_0 = Eigen::VectorXd::Zero(R), all_l2 = Eigen::VectorXd::Constant(max_models, 1e-6);
    Eigen::VectorXi freq_ind = Eigen::VectorXi::Zero(months);

    // Fits
    CrossVal::CV_res results(R);
    results.min_mse = DBL_MAX;

    // Create a vector indicating which variables are of the lowest frequency (aon only monthly)
    for (int n = 0; n < N; ++n)
    {
        if (frequency(n) == 12)
        {
            freq_ind(ind) = n;
            ++ind;
        }
    }

    // Initialize the random model sets
    std::uniform_real_distribution<> dis_l2(l2_min, l2_max);

    if (l1_min_nan || l1_max_nan) {

        /* Maximum and minim l1 penalty not provided -> Draw from the vector of all empirically relevant l1 penalties */

        // Decompositions
        JacobiSVD<Eigen::MatrixXd> SVD_i, X_SVD(X, ComputeThinU | ComputeThinV);
        ColPivHouseholderQR<Eigen::MatrixXd> QR;

        // Vectors
        Eigen::VectorXd y_temp = Eigen::VectorXd::Zero(T), all_l1_thresholds = Eigen::VectorXd::Zero(0);

        // Matrices
        Eigen::MatrixXd V = X_SVD.matrixV(), Sigma = X_SVD.singularValues().asDiagonal();


        // Initialise A = V(,0:k-1)
        Eigen::MatrixXd A = V(all, seq(0, R - 1));
         
        for (int i = 1; i < max_models; ++i) {
            all_l2(i) = std::pow(10, dis_l2(gen));
            for (int k = 0; k < R; ++k) {

                // Create target vector
                y_temp = X * A.col(k);
                
                // Compute all empirically relevant l1 penalties
                Eigen::VectorXd Temp = LARS<true>(y_temp, X, all_l2(i));
                all_l1_thresholds.conservativeResize(all_l1_thresholds.size() + Temp.size());
                all_l1_thresholds.tail(Temp.size()) = Temp;

            }
            std::uniform_int_distribution<> dist_l1(0, static_cast<int>(all_l1_thresholds.size()) - 1);
            rnd_li_given_l2 = all_l1_thresholds(dist_l1(gen));
            all_l1.col(i).setConstant(rnd_li_given_l2);
        }

    }else{
        std::uniform_real_distribution<double> dis_l1(0.0, 1.0);
        for (int i = 1; i < max_models; ++i) {
            all_l2(i) = std::pow(10, dis_l2(gen));
            rnd_li_given_l2 = dis_l1(gen);
            all_l1.col(i).setConstant(rnd_li_given_l2 * l1_max);
        }
    }

    /* Validation loop */

    if (log) { std::cout << "\nCurrently cross-validating. Please stand by." << std::endl; }

    // Serialise Eigen operations within the parallelisation according to https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
    Eigen::setNbThreads(1);

    // Parallel looping
#pragma omp parallel for schedule(dynamic) shared(results)
    for (int i = 0; i < max_models; ++i) {

        // Set local thread variables
        const double l2 = all_l2(i);
        const Eigen::VectorXd curr_l1_local = all_l1.col(i);

        // Reset the predictor matrix in case the decorrelating transformation has been used in the last step (locally!)
        Eigen::MatrixXd X_local = X_in;

        // Calculate the SDFM
        Filtering::KFS_fit DFM;
        SparseDFM::SDFMKFS(DFM, X_local(all, freq_ind), delay(freq_ind), Eigen::VectorXi::Constant(R, INT_MAX), R, order, decorr_errors, crit, "LARS", l2, curr_l1_local, 0., NAN, NAN, 1, 1000, INT_MIN, comp_null,
            0, 0.0001, conv_threshold, 0);
        if (decorr_errors)
        {
            X_local(all, freq_ind) = (X_local(all, freq_ind) * DFM.C.transpose()).eval();

        }

        // Compute the BIC
        double BIC_curr = std::log((1. / (double(months) * double(T)) * (X_local(Eigen::seq(DFM.order - 1, last), freq_ind) - (DFM.Lambda_hat * DFM.F(all, Eigen::seq(0, last - 1))).transpose()).squaredNorm())) + double((DFM.Lambda_hat.array() != 0).count()) * (std::log(double(months) * double(T))) / (double(months) * double(T));


        // Compare results between threads in critical environment
#pragma omp critical
        {
            if (BIC_curr < results.min_mse)
            {
                results.l2 = l2;
                results.l1 = curr_l1_local;
                results.min_mse = BIC_curr;
            }
        }
    }

    /* End validaiton loop */

    // Enable Eigen's internal parallelisation again
    Eigen::setNbThreads(0);

    // End of the loop
    if (log) { std::cout << "\rDone. The following parametrization is chosen: l2 = " << results.l2 << "; Number of zero entries for each factor: " << results.sel << "; min{BIC} = " << results.min_mse << '\n'; }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Print the duration if wished
    if (timer)
    {
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;
    }

    return results;
}
