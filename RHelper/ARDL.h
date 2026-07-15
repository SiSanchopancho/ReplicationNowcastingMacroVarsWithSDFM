/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright (C) 2024 Domenic Franjic
 *
 * This file is part of TwoStepSDFM.
 *
 * TwoStepSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * TwoStepSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.
 */


#pragma once

#ifndef ARDL
#define ARDL

 // Including external libraries
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <stdlib.h>
#include <math.h>
#include <string>


Rcpp::List runARDL(
	Rcpp::NumericVector target_variable,
	Rcpp::NumericVector predictor_variable,
	const unsigned max_target_lags,
	const unsigned max_predictor_lags,
	const std::string crit
);

#endif /* defined(ARDL) */