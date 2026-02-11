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

#ifndef DATA_HANDLING
#define DATA_HANDLING

#ifdef _MSC_VER // Check if the compiler is MSVC
#pragma warning(disable : 4996) // Disable MSVC warning due to strcpy
#endif

// Including external libraries
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#define mkdir _mkdir // Define mkdir to use _mkdir on Windows
#else
#include <unistd.h> // For access and mkdir on POSIX systems
#endif

// Including internal libraries
#include "Developer.h"

namespace DataHandle {

    void demean(Eigen::MatrixXd& X);
    void demean(Eigen::VectorXd& v);
    double var(const VectorXd& v_in);
    Eigen::MatrixXd cov(const Eigen::MatrixXd& X_in);
    Eigen::MatrixXd corr(const Eigen::MatrixXd& X_in);
    double corr(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
    double corr(const Eigen::Block<Eigen::MatrixXd, -1, 1, true>& x, const Eigen::Block<Eigen::MatrixXd, -1, 1, true>& y);
    void standardise(Eigen::MatrixXd& X);
    void standardise(Eigen::VectorXd& v);
    void standardise(Block<MatrixXd, -1, -1, false>& v_in);
    Eigen::MatrixXd kroneckerProd(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y);
    void removeRow(Eigen::MatrixXd& X, const int& t, const bool& conservative = true);
    void removeCol(Eigen::MatrixXd& X, const int& n, const bool& conservative = true);
    void removeElement(Eigen::VectorXd& x, const int& t, const bool& conservative);
    void dataSaveCSV(const Eigen::MatrixXd& X, const std::string& name, const std::string& type);
    Eigen::MatrixXd dataLoadCSV(const std::string& name);
    Eigen::MatrixXd dataLoadCSV2(const std::string& name);
    Eigen::MatrixXd dataLoadCSVDate(const std::string& name);
    std::string namer(
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
    );
    bool folderExists(const std::string& folderPath);
    bool createFolder(const std::string& folderPath);
    bool ensureFolderExists(const std::string& folderName);
    void displayLoadingBar(int progress);

};
#endif /* defined(DATA_HANDLING) */
