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

#ifdef _MSC_VER

#include <complex.h>

#define LAPACK_COMPLEX_CUSTOM

typedef _Dcomplex lapack_complex_double;
typedef _Fcomplex lapack_complex_float;

#define lapack_complex_float_real(z)       (real(z))
#define lapack_complex_float_imag(z)       (imag(z))

#endif // _MSC_VER

