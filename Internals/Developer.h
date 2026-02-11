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

#ifndef DEVELOPER
#define DEVELOPER

 // Including external libraries
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <Eigen>

// Including internal libraries

template<typename Derived>
void MPrint(const Eigen::MatrixBase<Derived>& X, const char* name = "Matrix", int t1 = 9, int n1 = 9, int t0 = 0, int n0 = 0)
{
  if (t1 > static_cast<int>(X.rows()) - 1)
  {
    t1 = static_cast<int>(X.rows()) - 1;
  }
  if (n1 > static_cast<int>(X.cols()) - 1)
  {
    n1 = static_cast<int>(X.cols()) - 1;
  }
  if (t0 > static_cast<int>(X.rows()) - 1)
  {
    t1 = 0;

    std::cout << "Starting index for t0 is bigger then Matrix size and has been set to zero!" << '\n';

  }
  if (n1 > static_cast<int>(X.cols()) - 1)
  {
    n0 = 0;

    std::cout << "Starting index for n0 is bigger then Matrix size and has been set to zero!" << '\n';

  }
  std::cout << '\n' << name << "(" << t0 << ':' << t1 << ',' << n0 << ':' << n1 << ") = \n" << X.block(t0, n0, t1 - t0 + 1, n1 - n0 + 1) << '\n';
}

template<typename Derived>
void Dim(const Eigen::MatrixBase<Derived>& X, const char* matrix = "Matrix")
{
  std::cout << '\n' << matrix << " is (" << X.rows() << "x" << X.cols() << ")" << '\n';
}

using namespace Eigen;

namespace Developer {

  extern

    // Development helpers
    void Print(const double&, const char*); // simple printing of double
    void iPrint(const int&, const char*); // simple printing of integer
    void Equal(const MatrixXd&, const MatrixXd&, const char* X_name = "X", const char* Y_name = "Y"); // prints whether RHS and LHS matrices are (approximately) equal
    void Equal(const VectorXd&, const VectorXd&, const char* x_name = "X", const char* y_name = "Y"); // prints whether RHS and LHS vector are (approximately) equal
    void here(); // error ind
    void here(const int&); // error ind for use within in loops
};
#endif /* defined(DEVELOPER) */
