#include "Orders.h"


double Orders::NoOfDF(
    const MatrixXd& X_in, // Data Matrix
    const int& k0, // lower bound for the factor estimates
    const int& k1, // upper bound for the factor estimates
    const double& w0, // Test frequency omega
    const VectorXi& s // Multiplicators to calculate the approximating frequencies
)
{
    // Estimating the number of factors according to Onatskiy (2009)

    // Dummies

    // Matrices
    MatrixXd X = X_in;

    // Integers
    int T = X.rows(), M = s.size();

    // Reals
    double R = DBL_MIN;

    // Vectors
    VectorXd v = VectorXd::LinSpaced(T, 1, T);

    // complex
    std::complex<double> i(0, 1);

    // Misshandling
    if (7 < k1 - k0) {
        std::cout << '\n' << "Warning! Power of the test might be low since k1 - k0 > 7." << '\n';
    }
    else if (18 < k1 - k0)
    {
        std::cout << '\n' << "Error! Critical values for k1 - k0 > 18 not available. Decrease k1." << '\n';
        return EXIT_FAILURE;
    }
    if (double(T) / 2. < M)
    {
        std::cout << '\n' << "Error! M must be strictly smaller then T/2" << '\n';
        return EXIT_FAILURE;
    }
    if (M < 5 * (k1 - k0))
    {
        std::cout << '\n' << "Warning! Power of the test might be low since M < 5 * (k1 - k0). Increase M and/or decrease k1." << '\n';
    }
    for (int m = 0; m < M; ++m)
    {
        for (int mm = m; mm < M; ++mm)
        {
            if ((m == mm) && (s(m) == 0 || s(m) == T / 2))
            {
                std::cout << '\n' << "Error! s_i cannot be 0 or T/2 for all s_i in s" << '\n';
                return EXIT_FAILURE;
            }
            else if ((s(m) + s(mm)) % T == 0)
            {
                std::cout << '\n' << "Error! (s(" << m << ") + s(" << mm << ")) % T == 0. This is not allowed for any two s_i and s_j in s." << '\n';
                return EXIT_FAILURE;
            }
        }
    }

    // Calculate descrete fourie transformation

    MatrixXcd F = (-i * 2. * M_PI * double(1. / T) * (v * s.transpose().cast<double>()).cast<std::complex<double>>()).array().exp().matrix();
    MatrixXcd Xf = F.transpose() * X;

    // Eigen decompositions
    MatrixXcd Gram = Xf * Xf.adjoint();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(Gram);
    VectorXd e = eig.eigenvalues().reverse();

    // Calculate the test statistic
    R = ((e(seq(k0, k1 - 1)) - e(seq(k0 + 1, k1))).array() / (e(seq(k0 + 1, k1)) - e(seq(k0 + 2, k1 + 1))).array()).maxCoeff();

    MatrixXd test = DataLoadCSV("./Imports/Onatski_test_stats_csv.txt");

    return double((test(all, k1 - k0 - 1).array() > R).count()) / 1000;
}

double Orders::NoOfSF(
    const MatrixXd& X_in, // Data Matrix
    const int& k0, // Lower bound for the factor estimates
    const int& k1 // Upper bound for the factor estimates
)
{
    // Estimating the number of factors according to Onatskiy (2009)

    // Dummies

    // Matrices
    MatrixXd X = X_in;

    // Integers
    int T = X.rows(), cutoff = std::floor(double(T) / 2.);

    // Reals
    double R = DBL_MIN;

    // Vectors
    VectorXd v = VectorXd::LinSpaced(T, 1, T);

    // complex
    std::complex<double> i(0, 1);

    // Misshandling
    if (7 < k1 - k0) {
        std::cout << '\n' << "Warning! Power of the test might be low since k1 - k0 > 7." << '\n';
    }
    else if (18 < k1 - k0)
    {
        std::cout << '\n' << "Error! Critical values for k1 - k0 > 18 not available. Decrease k1." << '\n';
        return EXIT_FAILURE;
    }

    // Transform the data
    MatrixXcd Xf = X(seq(0, cutoff - 1), all) + i * X(seq(cutoff, 2 * cutoff - 1), all);

    // Eigen decompositions
    MatrixXcd Gram = Xf * Xf.adjoint();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(Gram);
    VectorXd e = eig.eigenvalues().reverse();

    // Calculate the test statistic
    R = ((e(seq(k0, k1 - 1)) - e(seq(k0 + 1, k1))).array() / (e(seq(k0 + 1, k1)) - e(seq(k0 + 2, k1 + 1))).array()).maxCoeff();

    MatrixXd test = DataLoadCSV("./Imports/Onatski_test_stats_csv.txt");

    return double((test(all, k1 - k0 - 1).array() > R).count()) / 1000;
}

int Orders::NoOfFactors(
    const MatrixXd& X, // As in NoOFDF
    const int& k0, // As in NoOFDF
    const int& k1, // As in NoOFDF
    const double& w0, // As in NoOFDF
    const VectorXi& s, // As in NoOFDF
    const double& alpha, // Confidence level
    const bool& log // Talk to me
)
{
    // Estimating the number of factors according to Onatskiy (2009)

    // Integers
    int k = k0 - 1;

    // Reals
    double p = DBL_MIN;

    while (p < alpha && k < k1)
    {
        ++k;
        p = NoOfDF(X, k, k1, w0, s);
    }

    if (k + 1 == k1)
    {
        std::cout << '\n' << "Warning! k has been chosen as k1 - 1. It might be necessary to increase k1 and repeat the procedure" << '\n';
    }

    if (log) std::cout << '\n' << "k = " << k << " with a p-value of " << p << " and a critical value of alpha = " << alpha << '\n';

    return k;
}

int Orders::NoOfFactors(
    const MatrixXd& X, // As in NoOFDF
    const int& k0, // As in NoOFDF
    const int& k1, // As in NoOFDF
    const double& alpha,  // Confidence level
    const bool& log // Talk to me
)
{
    // Estimating the number of factors according to Onatskiy (2009)

    // Integers
    int k = k0 - 1;

    // Reals
    double p = DBL_MIN;

    while (p < alpha && k < k1)
    {
        ++k;
        p = NoOfSF(X, k, k1);
    }

    if (k + 1 == k1)
    {
        std::cout << '\n' << "Warning! k has been chosen as k1 - 1. It might be necessary to increase k1 and repeat the procedure" << '\n';
    }

    if (log) std::cout << '\n' << "k = " << k << " with a p-value of " << p << " and a critical value of alpha = " << alpha << '\n';

    return k;
}

int Orders::VARorder(
    const MatrixXd& F, // Data matrix
    const int& O, // Maximum order to be tested
    const double& comp_null, // Comptational zero
    const char* crit // Information Criterion used to determine the order ("BIC", "AIC", "HIC")
)
{
    // Infer the order of a VAR process using IC criteria

    // Dummies

    // Integers
    int K = F.rows(), T = F.cols(), order = 0;

    // Matrices
    MatrixXd curr = MatrixXd::Constant(O, 4, DBL_MAX);

    for (int o = 1; o <= O; ++o)
    {
        MatrixXd F_curr(K * o, T);

        // Collect all lags inside one matrix

        for (int oo = 0; oo < o; ++oo)
        {

            F_curr.block(K * oo, oo, K, T - oo) = F.block(0, 0, K, T - oo);

        }

        // OLS

        MatrixXd F_t = F_curr(all, seq(o + 1, last)).transpose();
        MatrixXd F_t_lag = F_curr(all, seq(o, last - 1)).transpose();
        MatrixXd Phi = ((F_t_lag.transpose() * F_t_lag).llt().solve(MatrixXd::Identity(K * o, K * o)) * F_t_lag.transpose() * F_t).transpose();
        Phi = Phi.unaryExpr([comp_null](double x) {return (comp_null < std::abs(x)) ? x : 0.; });
        MatrixXd F_hat = (Phi.topLeftCorner(K, K * o) * F_curr(all, seq(o, last - 1))).transpose();

        // Calculate the residuals variance-covariance-matrix and its determinant

        MatrixXd Res_Var_Cov = (F_t(0, seq(0, K - 1)) - F_hat.row(0)).transpose() * (F_t(0, seq(0, K - 1)) - F_hat.row(0));
        for (int t = 1; t < T - 1 - o; ++t)
        {
            Res_Var_Cov += (F_t(t, seq(0, K - 1)) - F_hat.row(t)).transpose() * (F_t(t, seq(0, K - 1)) - F_hat.row(t));
        }
        Res_Var_Cov /= double(T);

        double det = Res_Var_Cov.determinant();

        // Calculate the information criterion

        curr(o - 1, 0) = o;
        curr(o - 1, 1) = std::log(std::abs(det)) + 2. * double(K * K * o) / double(T); // AIC
        curr(o - 1, 2) = std::log(std::abs(det)) + double(K * K * o) * std::log(double(T)) / double(T); // BIC
        curr(o - 1, 3) = std::log(std::abs(det)) + 2. * double(K * K * o) * std::log(std::log(double(T))) / double(T); // HIC
    }

    //MPrint(curr, "curr");

    if (!strcmp(crit, "AIC"))
    {
        curr.col(1).minCoeff(&order);
        return order + 1;
    }
    else if (!strcmp(crit, "BIC"))
    {
        curr.col(2).minCoeff(&order);
        return order + 1;
    }
    else if (!strcmp(crit, "HIC"))
    {
        curr.col(3).minCoeff(&order);
        return order + 1;
    }
    else
    {
        std::cout << '\n' << "Error! Crit must be one of the following: \"BIC\", \"AIC\", or \"HIC\"." << '\n';
        return EXIT_FAILURE;
    }
}

double Orders::DMtest(
    const VectorXd& fe1, // Forecasting error of the first model
    const VectorXd& fe2 // Forecasting error of the second model
)
{

    // set up the test statistic

    const double T = double(fe1.size());
    VectorXd d = (fe1.array() * fe1.array()) - (fe2.array() * fe2.array());
    const double d_bar = d.mean();
    d -= d_bar * VectorXd::Ones(T);
    const double gamma_0 = (1. / T) * d.dot(d);

    // Calculate the sum of autocorrelations

    double gamma_k_sum = 0.;
    double gamma_k = 0.;
    const int H = std::ceil(std::pow(T, (1. / 3.)) + 1);
    for (int k = 0; k < H; ++k)
    {
        gamma_k = 0.;
        for (int i = k; i < T; ++i)
        {
            gamma_k += d(i) * d(i - k);
        }
        gamma_k *= (1. / T);
        gamma_k_sum += gamma_k;
    }
    gamma_k_sum *= 2;

    // Calculate the test-statistic

    const double DM = d_bar / (std::sqrt((gamma_0 + gamma_k_sum) / T));

    // Calculate the corresponding alpha value

    const double mean = 0.0;
    const double stddev = 1.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(mean, stddev);

    double alpha = 0.5 * (1.0 + std::erf((DM - mean) / (stddev * std::sqrt(2.0))));

    return alpha;

}