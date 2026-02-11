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

#include "DataHandle.h"

/* In-Place column de-meaning of a matrix */
void DataHandle::demean(Eigen::MatrixXd& X)
{
    Eigen::VectorXd mu = X.colwise().mean();
    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n);
        }
    }
    return;
};

/* In-Place de-meaning of a vecctor*/
void DataHandle::demean(Eigen::VectorXd& v)
{
    v -= v.mean() * Eigen::VectorXd::Ones(v.size());
    return;
};

/* Variance of a vector*/
double DataHandle::var(const VectorXd& v_in)
{
    double mean = v_in.mean();
    Eigen::VectorXd centered = v_in.array() - mean;
    return (centered.squaredNorm()) / (v_in.size() - 1);
};

/* Compute the variance co-variance matrix */
Eigen::MatrixXd DataHandle::cov(const Eigen::MatrixXd& X_in)
{
    Eigen::MatrixXd X = X_in;
    Eigen::VectorXd mu = X.colwise().mean();
    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n);
        }
    }
    return (X.transpose() * X) / double(X.rows() - 1);
};

/* Compute the correlation matrix */
Eigen::MatrixXd DataHandle::corr(const Eigen::MatrixXd& X_in)
{
    Eigen::MatrixXd X = X_in;
    Eigen::VectorXd mu = X.colwise().mean();
    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n);
        }
    }
    Eigen::MatrixXd cov = (X.transpose() * X) / double(X.rows() - 1);
    Eigen::VectorXd ss = cov.diagonal();
    Eigen::MatrixXd S_sqrt_inv = Eigen::MatrixXd::Zero(X.cols(), X.cols());
    for (int n = 0; n < X.cols(); ++n)
    {
        S_sqrt_inv(n, n) = 1 / sqrt(ss(n));
    }
    return S_sqrt_inv * cov * S_sqrt_inv;
};

/* Compute the correlation between two vectors */
double DataHandle::corr(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
    Eigen::VectorXd x_in = x;
    Eigen::VectorXd y_in = y;
    demean(x_in);
    demean(y_in);
    return (x_in.dot(y_in)) / (x_in.norm() * y_in.norm());
};

/* As above but for Eigen::Block */
double DataHandle::corr(const Eigen::Block<Eigen::MatrixXd, -1, 1, true>& x, const Eigen::Block<Eigen::MatrixXd, -1, 1, true>& y)
{
    Eigen::VectorXd x_in = x;
    Eigen::VectorXd y_in = y;
    demean(x_in);
    demean(y_in);
    return (x_in.dot(y_in)) / (x_in.norm() * y_in.norm());
};

/* In-place column-wise data standardisation*/
void DataHandle::standardise(Eigen::MatrixXd& X)
{
    Eigen::VectorXd mu = X.colwise().mean();
    Eigen::VectorXd sigma = cov(X).diagonal();
    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n);
            X(t, n) /= std::sqrt(sigma(n));
        }
    }
    return;
};

/* In-Place vector standardisation*/
void DataHandle::standardise(Eigen::VectorXd& v)
{
    v -= v.mean() * Eigen::VectorXd::Ones(v.size());
    double s = std::sqrt(var(v));
    for (int n = 0; n < v.size(); ++n)
    {
        v(n) /= s;
    }
    return;
};

/* As above but for Eigen::Block */
void DataHandle::standardise(Block<MatrixXd, -1, -1, false>& v_in)
{
    VectorXd v = v_in;
    v -= v.mean() * VectorXd::Ones(v.size());
    double s = std::sqrt(var(v));
    for (int n = 0; n < v.size(); ++n)
    {
        v(n) /= s;
    }
    v_in = v;
    return;
};

/* Kronecker product*/
Eigen::MatrixXd DataHandle::kroneckerProd(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y)
{
    const int M = static_cast<int>(X.rows());
    const int N = static_cast<int>(X.cols());
    const int P = static_cast<int>(Y.rows());
    const int R = static_cast<int>(Y.cols());
    Eigen::MatrixXd K(M * P, N * R);
    for (int n = 0; n < N; ++n)
    {
        for (int m = 0; m < M; ++m)
        {
            K.block(n * P, m * R, P, R) = X(m, n) * Y;
        }
    }
    return K;
};


/* Remove row from matrix in-place */
void DataHandle::removeRow(Eigen::MatrixXd& X, const int& t, const bool& conservative)
{
    int T = static_cast<int>(X.rows()) - 1;
    int N = static_cast<int>(X.cols());
    if (t < T)
    {
        X.block(t, 0, T - t, N) = X.block(t + 1, 0, T - t, N).eval();
    }
    if (conservative)
    {
        X.row(T).setZero();
    }
    else
    {
        X.conservativeResize(T, N);
    }

};

/* Remove row from matrix in-place*/
void DataHandle::removeCol(Eigen::MatrixXd& X, const int& n, const bool& conservative)
{
    int T = static_cast<int>(X.rows());
    int N = static_cast<int>(X.cols()) - 1;
    if (n < N)
    {
        X.block(0, n, T, N - n) = X.block(0, n + 1, T, N - n).eval();
    }
    if (conservative)
    {
        X.col(N).setZero();
    }
    else
    {
        X.conservativeResize(T, N);
    }
};

/* Remove element from vector in-place */
void DataHandle::removeElement(Eigen::VectorXd& x, const int& t, const bool& conservative)
{
    int T = static_cast<int>(x.size()) - 1;

    if (t < T)
    {
        x.block(t, 0, T - t, 1) = x.block(t + 1, 0, T - t, 1).eval();
    }

    if (conservative)
    {
        x(T) = 0.;
    }
    else
    {
        x.conservativeResize(T, 1);
    }
};

/* Save matrix as csv */
void DataHandle::dataSaveCSV(const Eigen::MatrixXd& X, const std::string& name, const std::string& type)
{
    std::string out_name = name + type;
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file;
    file.open(out_name, std::ofstream::out | std::ofstream::trunc);
    if (file.is_open())
    {
        file << X.format(CSVFormat);
    }
    file.close();
};

/* Load .csv to Eigen::MatrixXd */
Eigen::MatrixXd DataHandle::dataLoadCSV(const std::string& name)
{
    std::vector<double> entry;
    std::string rows, element;
    std::ifstream matrixDataFile(name);
    int t = 0;
    while (std::getline(matrixDataFile, rows))
    {
        std::stringstream matrixRowStringStream(rows);

        while (std::getline(matrixRowStringStream, element, ','))
        {
            entry.push_back(std::stod(element));
        }
        ++t;
    }
    return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(entry.data(), t, entry.size() / t);
};

/* Load .csv to Eigen::MatrixXd ColMajor-wise */
Eigen::MatrixXd DataHandle::dataLoadCSV2(const std::string& name)
{
    std::vector<double> entry;
    std::string rows, element;
    std::ifstream matrixDataFile(name);
    int t = 0;
    while (std::getline(matrixDataFile, rows))
    {
        std::stringstream matrixRowStringStream(rows);

        while (std::getline(matrixRowStringStream, element, ','))
        {
            entry.push_back(std::stod(element));
        }
        ++t;
    }
    return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(entry.data(), t, entry.size() / t);
};

/* Load in "dates" */
MatrixXd DataHandle::dataLoadCSVDate(const std::string& name)
{
    std::vector<std::string> dateStrings;
    std::string rows;
    std::ifstream matrixDataFile(name);
    while (getline(matrixDataFile, rows))
    {
        dateStrings.push_back(rows);
    }
    std::vector<std::tm> dates;
    for (const auto& dateString : dateStrings)
    {
        std::tm date = {};
        std::istringstream ss(dateString);
        ss >> std::get_time(&date, "%Y-%m-%d");
        if (ss.fail()) {
            std::cerr << "Error parsing date: " << dateString << std::endl;
        }
        dates.push_back(date);
    }
    Eigen::MatrixXd result(dates.size(), 3);
    for (size_t i = 0; i < dates.size(); ++i)
    {
        result(i, 0) = dates[i].tm_year + 1900;  // Years since 1900
        result(i, 1) = dates[i].tm_mon + 1;      // Months are 0-based
        result(i, 2) = dates[i].tm_mday;    // Months are 0-based
    }
    return result;
}

/* Convenient namer for storing the results */
std::string DataHandle::namer(
    const std::string& folder_name,
    const int& K_sim,
    const int& N,
    const int& T,
    const double& prob,
    const double& beta_param,
    const bool& corr,
    const int FCH,
    const int& T_CV,
    const int& G_rnd,
    const int& seed
)
{
    const size_t size = folder_name.size() + 100;
    std::string out_name;
    out_name.reserve(size);
    out_name = folder_name + "_K" + std::to_string(K_sim) +
        "_N" + std::to_string(N) +
        "_T" + std::to_string(T) +
        "_p" + std::to_string(static_cast<int>(prob * 100)) +
        "_b" + std::to_string(static_cast<int>(std::round(beta_param * 1000))) +
        "_c" + std::to_string(corr) +
        "_H" + std::to_string(FCH) +
        "_TCV" + std::to_string(T_CV) +
        "_G1" + std::to_string(G_rnd) +
        "_seed" + std::to_string(seed);
    return out_name;
};

/* Check whether the folder exists for MSVC (top) and g++ (bottom) */
bool DataHandle::folderExists(const std::string& folderPath) {
#ifdef _WIN32
    struct _stat info;
    return _stat(folderPath.c_str(), &info) == 0 && (info.st_mode & _S_IFDIR);
#else
    struct stat info;
    return stat(folderPath.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
#endif
};

/* Create folder conditionally */
bool DataHandle::createFolder(const std::string& folderPath) {
#ifdef _WIN32
    return _mkdir(folderPath.c_str()) == 0;
#else
    return mkdir(folderPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0;
#endif
};

/* Wrapper for the above functions */
bool DataHandle::ensureFolderExists(const std::string& folderName) {
    if (!folderExists(folderName)) {
        if (createFolder(folderName)) {
            std::cout << "Folder '" << folderName << "' created successfully." << std::endl;
            return true;
        }
        else {
            std::cerr << "Error creating folder '" << folderName << "'." << std::endl;
            return false;
        }
    }
    else {
        std::cout << "Folder '" << folderName << "' already exists." << std::endl;
        return true;
    }
};

/* Does what it says on the package*/
void DataHandle::displayLoadingBar(int progress) {

    const char chars[] = { '-', '/', '|', '\\' };
    const int numChars = sizeof(chars) / sizeof(chars[0]);

    std::cout << "\rProgress: " << chars[progress % numChars] << " " << std::flush;
}