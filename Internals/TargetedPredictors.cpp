#include "TargetedPredictors.h"

CrossVal::CV_res TargetedPredictors::TPCrossValDFM(
    const int& x0_ind, // Index of the variable of interest
    const MatrixXd& X, // (T x N) Data matrix
    const int& H, // Number of forecasts for the cross-validation
    const int& K, // Number of factors of the DFm
    const VectorXi& sel, // Hyperparameter grid (as of now n umber of selected variables)
    const VectorXd& l2, // Vector of potential l2 penalties
    const VectorXi& delay, // Vector of delays
    const VectorXi& date, // date vector
    const VectorXi& frequency, // frequency vector
    const bool& decorr_errors, // Decorrelate the idiosyncratic errors
    const int& order, // Max order of the VAR process to be checked in the order estimation
    const char* crit, // Cross-validation method
    const double& comp_null, // Computational zero
    const double& conv_threshold, // Conversion threshold for the SPCA algorithm
    const bool& timer
)
{

    // Cross-Validation Wrappper

    // Start timer

    auto start = high_resolution_clock::now();

    // Dummies

    // Integers
    int ind = 0, T = X.rows(), N = X.cols(), len_sel = sel.rows(), months = (frequency.array() == 12).count(), goal = std::pow(len_sel, K), running = 0;

    // Reals
    double progress = 0, progr = 0.;

    // Matrices
    MatrixXd IdentK = MatrixXd::Identity(K, K);

    // Vectors
    VectorXd lambda_0 = VectorXd::Zero(K);
    VectorXi ind_vec = VectorXi::Zero(K), freq_ind = VectorXi::Zero(months), fc_dates = VectorXi::Zero(H);

    // Fits
    CrossVal::CV_res results(1);
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

        fc_dates.resize(std::floor(double(H) / 3.));

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

        fc_dates = VectorXi::LinSpaced(H, T - H, T - 1);

    }

    // Print the fc dates for monitoring purposes

    int cv_dates = fc_dates.size();

    MatrixXi Check_dates = MatrixXi(int(floor(double(cv_dates) / 2.)), 2);
    Check_dates.col(0) = fc_dates.head(int(floor(double(cv_dates) / 2.)));
    Check_dates.col(1) = fc_dates.tail(int(floor(double(cv_dates) / 2.)));

    //MPrint(Check_dates, "CV_Window_Start_End");

    // Create dummies to store the fes and squares of the residuals

    MatrixXd fes = MatrixXd::Zero(fc_dates.size(), 2);
    VectorXd res_sq = VectorXd::Zero(fc_dates.size());

    // Main cv loop

    // The outer loop loops over all l2 values while the inner loop loops over all l1 combinations given the number of factors k
    // Overall, the CV loop will check |l2| * |l1| ^ K (aon)

    for (int i = 0; i < l2.size(); ++i)
    {
        // Give some feedback on what your doing

        std::cout << '\n' << "Cross validating combination for sel given l2 = " << l2(i) << '\n' << "Please stand by." << "      ";

        // Set the vector indicating the current l1 parameter combination to 0

        ind_vec.setZero();

        progr = 0.;

        for (int curr_sel : sel)
        {

            // Calculate the state-space model for the current hyperparameters

            // Internal constants

            const int max_iterations = 1000;

            // Progressbar

            progress = std::round(100. * double(progr) / double(goal));

            if (std::floor(double(progress) / 10.) == 0.)
            {
                std::cout << "\b\b\b" << progress << "% " << std::flush;
            }
            else if (std::floor(double(progress - 1) / 10.) == 0. && 1. <= double(progress) / 10.)
            {
                std::cout << "\b\b\b" << progress << "% " << std::flush;
            }
            else
            {
                std::cout << "\b\b\b\b" << progress << "% " << std::flush;
            }

            // Calculate the oos forecasting errors

            int ind2 = 0;
            running = 0;
            for (int h : fc_dates)
            {

                // Progress animation

                if (running % 4 == 0)
                {
                    std::cout << " -" << std::flush;
                }
                else if (running % 4 == 1)
                {
                    std::cout << " \\" << std::flush;
                }
                else if (running % 4 == 2)
                {
                    std::cout << " |" << std::flush;
                }
                else if (running % 4 == 3)
                {
                    std::cout << " /" << std::flush;
                }
                ++running;

                // Create a data matrix of aggregated monthly variables

                VectorXi sel_ind;
                VectorXi del_TP;

                TargetedPredictors::Selector(sel_ind, del_TP, X, delay, frequency, h, x0_ind, 5, NAN, curr_sel);

                KFS_fit DFM;
                DFMKFS(DFM, X(seq(0, h), sel_ind), del_TP, K, 10, decorr_errors);

                // Nowcast the variable of interest

                // Calculate the OLS fit of the loadings of the VOI on the estimated factors

                // Calculate and save the nowcast

                fes(ind2, 1) = Forecast::factorForecaster(X, DFM, K, frequency(x0_ind), delay(x0_ind), h, x0_ind);
                ++ind2;
                std::cout << "\b\b" << std::flush;
            }

            // Save the residuals squares

            res_sq = (fes.col(1).array() * fes.col(1).array()).matrix();

            // Save the results of this run if it is the mean of the squared residuals is smaller then the last best run. Else, discard the results.

            if (res_sq.mean() < results.min_mse)
            {
                results.l2 = l2(i);
                results.sel(0) = curr_sel;
                res_sq = (fes.col(1).array() * fes.col(1).array()).matrix();
                results.min_mse = res_sq.mean();
                fes.col(0) = fes.col(1);
            }

            // Go to the next K-ple of hyper-parameters

            ++ind_vec(K - 1);

            if (len_sel - 1 < ind_vec(K - 1) && 1 < K)
            {
                for (int k = K - 1; 0 < k; --k)
                {
                    if (len_sel - 1 < ind_vec(k))
                    {
                        ind_vec(k) = 0;
                        ++ind_vec(k - 1);
                    }
                }
            }

            ++progr;
        }

        // Reset progressbar
        std::cout << "\b\b\b\b\b" << std::flush;
        std::cout << "Done." << '\n' << '\n';

        std::cout << "Targeted Predictors Crossvalidation completed. The current minimum CV-error is " << results.min_mse << " for l2 = " << results.l2 << " and sel = " << results.sel << "." << '\n' << '\n';
    }

    // End timer

    auto end = high_resolution_clock::now();

    // Calculate the duration in milliseconds

    auto duration = duration_cast<milliseconds>(end - start);

    if (timer)
    {

        // Print the duration
        std::cout << "TP CV took " << duration.count() << " milliseconds" << std::endl;

    }

    return results;
}

void TargetedPredictors::Selector(
    VectorXi& sel_ind,
    VectorXi& del_TP,
    const MatrixXd& X,
    const VectorXi& delay,
    const VectorXi& frequency,
    const int h,
    const int x0_ind,
    const double& l2,
    double l1,
    int selected_in,
    int steps,
    const double& comp_null,
    const bool& log
)
{
    // Quartify data matrix.

    // Dummies

    // Integers
    int N = X.cols(), T = X.rows();

    // Aggregate the observations
    int last_date = delay(x0_ind);
    MatrixXd Temp = X(seq(0, last - last_date), all);
    Temp.conservativeResize(X(seq(0, last - last_date), all).rows() + 2, X(seq(0, last - last_date), all).cols());
    Temp.setZero();
    Temp.bottomLeftCorner(T - last_date, N) = X(seq(0, last - last_date), all);

    MatrixXd Xq = MatrixXd::Zero((h - last_date + 1) / 3, N);
    for (int n = 0; n < N; ++n)
    {

        int t = 4;
        int tt = 0;

        if (frequency(n) == 4)
        {
            for (int ttt = 0; ttt <= h - last_date; ++ttt)
            {
                if ((ttt + 1) % 3 == 0)
                {
                    Xq(tt, n) = X(ttt, n);
                    ++tt;
                }
            }
        }
        else
        {
            while (t <= h - last_date + 2)
            {

                Xq(tt, n) = 1. / 3. * (Temp(t, n) + 2. * Temp(t - 1, n) + 3. * Temp(t - 2, n) + 2. * Temp(t - 3, n) + Temp(t - 4, n));
                t += 3;
                ++tt;
            }
        }
    }

    removeRow(Xq, 0, 0);

    // Estimate the model using only the monthly observaitons

    if (x0_ind != 0)
    {
        std::cout << '\n' << "Error! The variable of interest must be stored in the first column of X!" << '\n';

        return ;
    }
    
    MatrixXd Beta = LARS<false>(Xq.col(x0_ind), Xq(all, seq(1, last)), l2, l1, selected_in, steps, comp_null);
    
    sel_ind = VectorXi::Zero((Beta.array() != 0).count());
    del_TP = VectorXi::Zero((Beta.array() != 0).count());

    int ind = 0;
    for (int n = 1; n < N; ++n)
    {
        if (Beta(n - 1, 0) != 0)
        {
            sel_ind(ind) = n;
            del_TP(ind) = delay(n);
            ++ind;
        }
    }

    return;
}