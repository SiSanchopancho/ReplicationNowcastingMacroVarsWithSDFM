/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
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

#include "Developer.h"

////////////// Development Helper ///////////////////////////////////////////////////////////////////////////////////////////////////////////

void Developer::Print(const double& obj, const char* print)
{
	std::cout << '\n' << print << " =" << obj << '\n';
}

void Developer::iPrint(const int& obj, const char* print)
{
	std::cout << '\n' << print << " =" << obj << '\n';
}

void Developer::Equal(const MatrixXd& X, const MatrixXd& Y, const char* X_name, const char* Y_name)
{
	std::cout << '\n' << X_name << " == " << Y_name << " ?: " << X.isApprox(Y, 10e-10) << '\n';
}

void Developer::Equal(const VectorXd& x, const VectorXd& y, const char* x_name, const char* y_name)
{
	std::cout << '\n' << x_name << " == " << y_name << " ?: " << x.isApprox(y, 10e-10) << '\n';
}

void Developer::here()
{
    std::cout << '\n' << "There is a problem here." << '\n';
}

void Developer::here(const int& i)
{
    std::cout << '\n' << "There is a problem at iteration step " << i << "." << '\n';
}

