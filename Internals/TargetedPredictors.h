#pragma once

#ifndef TP
#define TP

// Including external libraries

#include "MSCVComplexHandler.h"
#include <lapacke.h>
#include <Eigen>
#include <stdlib.h>
#include <math.h>
#include <chrono>

// Including internal libraries

#include "Developer.h"
#include "Filtering.h"
#include "SparsePCA.h"
#include "SparseDFM.h"
#include "CrossVal.h"

using namespace Eigen;
using namespace Developer;
using namespace DataHandle;
using namespace SparsePCA;
using namespace Filtering;
using namespace SparseDFM;
using namespace CrossVal;
using namespace std::chrono;

namespace TargetedPredictors {

    extern
        CrossVal::CV_res TPCrossValDFM(
            const int& x0_ind,
            const MatrixXd& X,
            const int& H,
            const int& K,
            const VectorXi& sel,
            const VectorXd& l2,
            const VectorXi& delay,
            const VectorXi& date,
            const VectorXi& frequency,
            const bool& decorr_errors = 0,
            const int& order = 10,
            const char* crit = "BIC",
            const double& comp_null = 10e-15,
            const double& conv_threshold = 10e-6,
            const bool& timer = 1
        );

    void Selector(
        VectorXi& sel_ind,
        VectorXi& del_TP,
        const MatrixXd& X_in,
        const VectorXi& delay,
        const VectorXi& frequency,
        const int h,
        const int x0_ind,
        const double& l2 = 10e-6,
        double l1 = NAN,
        int selected_in = INT_MIN,
        int steps = INT_MIN,
        const double& comp_null = 10e-15,
        const bool& log = 0);

};
#endif /* defined(TP) */

