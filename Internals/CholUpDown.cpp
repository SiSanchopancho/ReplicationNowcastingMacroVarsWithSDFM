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

#include "CholUpDown.h"

/* Updating the Cholesky decomposition of the Gram matrix */
/* Source:
- Normal Splines. (2019, February 25). Algorithms for updating the Cholesky factorization. Normal Splines Blog. (https://normalsplines.blogspot.com/2019/02/algorithms-for-updating-cholesky.html)
*/
void CholUpDown::cholUpdate(
    Eigen::MatrixXd& L, // Lower triangular of the CD that is outh to be updated
    const Eigen::VectorXd& c, // Column of the data matrix as a vector that is added or removed to or from the data matrix
    const Eigen::MatrixXd& X, // Old data matrix (without the new column)
    const int& index, // Index of the column that has been removed from the old data matrix (-1 in case of upgrade)
    const double& l2 // l2 value of the elastic net problem
)
{

    /* Dummies */

    // Integers
    const int N = X.cols() + 1;

    //Vectors
    Eigen::VectorXd e = Eigen::VectorXd::Zero(N);

    //Matrices
    Eigen::MatrixXd XT_c = (1 / (1 + l2)) * (X.transpose() * c);

    // Calculated the updated lower triangular
    if (N == 2)
    {
        e(0) = XT_c(0, 0) / L(0, 0);
    }
    else
    {
        e.head(N - 1) = L.topLeftCorner(N - 1, N - 1).triangularView<Lower>().solve(XT_c);
    }

    e(N - 1) = sqrt((1 / (1 + l2) * (c.transpose() * c + l2)) - (e.head(N - 1).transpose() * e.head(N - 1))(0, 0));
    L.col(N - 1).head(N - 1) = Eigen::VectorXd::Zero(N - 1);
    L.row(N - 1).head(N) = e;

    return;
}

/* Downdating the Cholesky decomposition of the Gram matrix */
/* Source:
- Normal Splines. (2019, February 25). Algorithms for updating the Cholesky factorization. Normal Splines Blog. (https://normalsplines.blogspot.com/2019/02/algorithms-for-updating-cholesky.html)
*/
void CholUpDown::cholDowndate(
    Eigen::MatrixXd& L, // Lower triangular of the CD that is outh to be updated
    const Eigen::VectorXd& c, // Column of the data matrix as a vector that is added or removed to or from the data matrix
    const Eigen::MatrixXd& X, // Old data matrix (without the new column)
    const int& index, // Index of the column that has been removed from the old data matrix (-1 in case of upgrade)
    const double& l2 // l2 value of the elastic net problem
)
{

    /* Dummies */

    // Integers
    int NN = X.cols();

    // Reals
    double b = 1.0;

    // Vectors
    Eigen::VectorXd l = L(seq(index + 1, NN), index);

    //Remove the corresponding rows and columns of L
    DataHandle::removeRow(L, index);
    DataHandle::removeCol(L, index);
    Eigen::MatrixXd L_DL = L.block(index, index, NN - index, NN - index);
    int M = L_DL.cols();
    Eigen::VectorXd omega = l;

    // Update L recursively
    for (int m = 0; m < M; ++m)
    {

        double l_mm = std::sqrt(pow(L_DL(m, m), 2) + (1 / b) * pow(omega(m), 2));
        double gamma = pow(L_DL(m, m), 2) * b + pow(omega(m), 2);
        omega.tail(M - m - 1) -= (omega(m) / L_DL(m, m)) * L_DL.col(m).tail(M - m - 1);
        L(seq(index + m + 1, index + M - 1), index + m) = (l_mm / L_DL(m, m)) * L_DL(seq(m + 1, M - 1), m) + omega.tail(M - m - 1) * ((l_mm * omega(m)) / gamma);
        L(index + m, index + m) = l_mm;
        b += (pow(omega(m), 2) / pow(L_DL(m, m), 2));

    }

    return;
}