#include "ElNetSolve.h"

VectorXd ElNetSolve::LogSeq(const double& start, const double& end, const int& length)
{
	// Create logarithmically spaced sequence of lasso penalties

    // Dummies 

    // Vectors 

	VectorXd lseq = VectorXd::Zero(length);

    // Create the sequence

	double curr = start;
	double step = std::pow((end / start), double(1.0 / double(length - 1)));
	for (int i = 0; i < length; ++i)
	{
		lseq(i) = curr;
		curr *= step;
	}
	return lseq;
}

void ElNetSolve::CholUpdate(
    MatrixXd& L, // Lower triangular of the CD that is outh to be updated
    const VectorXd& c, // Column of the data matrix as a vector that is added or removed to or from the data matrix
    const MatrixXd& X, // Old data matrix (without the new column)
    const int& index, // Index of the column that has been removed from the old data matrix (-1 in case of upgrade)
    const double& l2, // l2 value of the elastic net problem
    const char* revision // Inicator for upgrade ("up") or downgrade ("down")
)
{

	// Updating/Downdating the Cholesky decomposition of the Gram matrix
	
    if (!strcmp(revision, "up"))
	{
		
        // Updating the Cholesky decomposition
		
        //Dummies

        // Integers

        const int N = X.cols() + 1;
		
        //Vectors

        VectorXd e = VectorXd::Zero(N);

        //Matrices

		MatrixXd XT_c = (1 / (1 + l2)) * (X.transpose() * c);

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
		L.col(N - 1).head(N - 1) = VectorXd::Zero(N - 1);
		L.row(N - 1).head(N) = e;

	}
	else if (!strcmp(revision, "down"))
	{

		// downdating the cholesky decomposition
		
        // Dummies

        // Integers

        int NN = X.cols();

        // Reals

        double b = 1;

        // Vectors

		VectorXd l = L(seq(index + 1, NN), index);
		
        //Remove the corresponding rows and columns of L

        removeRow(L, index);
		removeCol(L, index);
		MatrixXd L_DL = L.block(index, index, NN - index, NN - index);
		int M = L_DL.cols();
		VectorXd omega = l;
		
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
	}
}

MatrixXd ElNetSolve::LARS(
    const VectorXd& y,  // VOI
    const MatrixXd& X, // Predictors
    const double& l2, // l2 Penalty
    double l1, // l1 penalty (used for stopping)
    int selected_in, // Number of selected variables
    int steps, // Number of steps until stopping
    const double& comp_null, // Computational zero
    const bool& log // Talk to me
)
{

	// LARS Algorithmus (source elasticnet package R, spca package R, and original paper)

	// Dummies

	// Integers

    const int N = X.cols(), T = X.rows();
    int selected = selected_in, just_left = -1, curr_var_index = -1, gamma_index = -1, size_A_C = N,
        size_A = 0, s = 0, drop_index = -1;

    // Bools

	bool drop = 0;

    // Reals

	double gamma = DBL_MAX, gamma_hat = DBL_MAX, A_A = DBL_MAX, gamma_tilde = DBL_MAX, C_hat = DBL_MAX, l2_sqrt = sqrt(l2),
        l2_sqrt_inv = 1 / sqrt(1 + l2);

	// Matrices

	MatrixXd L = MatrixXd::Zero(N, N);

	// Vectors

	VectorXd w_A_vec = VectorXd::Zero(N), L_T_G_inv = VectorXd::Zero(N), G_A_inv_one = VectorXd::Zero(N), g_tilde = VectorXd::Zero(N),
        beta_curr = VectorXd::Zero(N), u_A = VectorXd::Zero(N + T), a(N), residuals = VectorXd::Zero(T + N),
        c = l2_sqrt_inv * X.transpose() * y, beta_hat = VectorXd::Zero(N), sign = VectorXd::Zero(N);
	VectorXi A_C_set = VectorXi::LinSpaced(N, 0, N - 1), A_set = -1 * VectorXi::Ones(N);

	// Special cases and misshandling

    if (l1 == 0. || l1 == INT_MIN)
    {
        l1 = NAN;
    }
	if (!(std::isnan(l1)) && (selected != INT_MIN))
	{
		std::cout << '\n' << "Error! Both penalty and number of variables is provided. Termination criterion is ambiguous. Only provide either of them." << '\n';
		return EXIT_FAILURE * VectorXd::Ones(N);
	}
	else if (l1 < 0 || (selected < 0 && selected != INT_MIN) || (steps < 0 && steps != INT_MIN))
	{
		std::cout << '\n' << "Error! A termination criterion ('l1', 'selected', or 'steps') has negative value." << '\n';
		return EXIT_FAILURE * VectorXd::Ones(N);
	}
	else if ((!(std::isnan(l1)) && ((2 * c.cwiseAbs().col(0).maxCoeff() / l2_sqrt_inv) <= l1)) || selected == 0)
	{
		std::cout << '\n' << "Warning! Either 'selected' is set to 0 or l1 is bigger or equal to smallest value for which the coefficient vector is the zero vector." << '\n';
		return VectorXd::Zero(N);
	}

	// Set the default stopping criterion with respect to number of variables selected
	
    if (l2 == 0.)
	{
		selected = (N <= (T - 1)) ? N : T - 1;
	}
	else if (N < selected)
	{
		selected = N;
	}

	// Redefinitions
	
    steps = ((steps == INT_MIN) ? 50 * (T <= N - 1 ? T : (N - 1)) : steps);
	residuals.head(T) = y;
	residuals.tail(N).setZero();
	VectorXd pen(steps + 1);
	pen(0) = c.cwiseAbs().col(0).maxCoeff();

	// Actual LARS loop
    
    while (s < steps && size_A < selected)
	{

		// Calculate current maximum correlation and retrieve the index of the corresponding variable evaluated at the current step
		
        C_hat = c(A_C_set.head(size_A_C)).cwiseAbs().maxCoeff(&just_left);

		if (log)  Print(c.cwiseAbs().col(0).sum(), "sum(|c|)");

		if (drop == 0)
		{

			// Previously no variable has just been dropped

			// Add the variable to the active set
			// Erase the variable from the inactive set
			
            curr_var_index = A_C_set(just_left);
			--size_A_C;
			++size_A;
			A_set(size_A - 1) = curr_var_index;
			sign(size_A - 1) = double(0 < c(curr_var_index)) - double(c(curr_var_index) < 0);
			A_C_set(seq(just_left, N - 2)) = A_C_set(seq(just_left + 1, N - 1));
			A_C_set(N - 1) = -1;

			if (s == 0)
			{

				// Create Gramm-Matrix
				
                L(0, 0) = sqrt(((X.col(curr_var_index).transpose() * X.col(curr_var_index))(0, 0) + l2) / (1 + l2));
			
            }
			else
			{

				// Update Cholesky
				
                CholUpdate(L, X(all, curr_var_index), X(all, A_set.head(size_A - 1)), -1, l2);
			
            }
		}
		else if (drop == 1 || size_A == N)
		{

			// A variable has just been dropped
            // Therefore, no new variable will be added since the last step has led to taking an "incomplete" LARS step
			
            curr_var_index = -1;
			drop = 0;
		
        }

		// Calculate equiengular vector
		
        G_A_inv_one.head(size_A) = sign.head(size_A);
		L.topLeftCorner(size_A, size_A).triangularView<Lower>().solveInPlace(G_A_inv_one.head(size_A));
		L.topLeftCorner(size_A, size_A).transpose().triangularView<Upper>().solveInPlace(G_A_inv_one.head(size_A));
		A_A = 1.0 / (std::sqrt((G_A_inv_one.head(size_A).transpose() * sign.head(size_A))(0)));
		VectorXd w_A = A_A * G_A_inv_one.head(size_A);
		w_A_vec(A_set.head(size_A)) = w_A;
		u_A.head(T) = l2_sqrt_inv * (X(all, A_set.head(size_A)) * w_A);
		u_A.tail(N) = l2_sqrt * l2_sqrt_inv * w_A_vec;
		
		// Calculate the maximum feasible step size (LARS-LASSO-Modification)
		
        gamma_tilde = DBL_MAX;

		for (int nn : A_set.head(size_A))
		{
			
            double g_curr = (-1 * (beta_hat(nn) / w_A_vec(nn)));
			
            if (g_curr < gamma_tilde && comp_null < g_curr)
			{

				gamma_tilde = g_curr;
				gamma_index = nn;
			
            }
		}

		if (N == size_A)
		{

            // For the last step, just go all the way and set gamma_hat to the maximum correlation

			gamma_hat = C_hat / A_A;
		
        }
		else
		{

			// Compute step size of the current LARS step
			
            a.head(size_A_C) = (X(all, A_C_set.head(size_A_C)).transpose() * u_A.head(T) + l2_sqrt * u_A.tail(N)(A_C_set.head(size_A_C))) * l2_sqrt_inv;
			gamma_hat = C_hat / A_A;

			for (int nn = 0; nn < size_A_C; ++nn)
			{

				double CAm = (C_hat - c(A_C_set(nn))) / (A_A - a(nn));
				double CAp = (C_hat + c(A_C_set(nn))) / (A_A + a(nn));

				if (CAm < gamma_hat && comp_null < CAm)
				{
					gamma_hat = CAm;
				}
				else if (CAp < gamma_hat && comp_null < CAp)
				{
					gamma_hat = CAp;
				}
			}
		}

		// Updating beta_hat, the residuals, and the correlation vector
        
		// If famma_hat is bigger then gamma_tilde, it is not possible to do a complete step due to sign restrictions on the coefficients
        // In this case the variable that would switch signs is dropped and the step is only as large as necessary to make the coefficient
        // of the variable that is dropped go to zero.

        gamma = (gamma_tilde < gamma_hat) ? gamma_tilde : gamma_hat;
		beta_curr = beta_hat;
		beta_hat += (gamma * w_A_vec);
		residuals -= (gamma * u_A);
		c = (X.transpose() * residuals.head(T) + l2_sqrt * residuals.tail(N)) * l2_sqrt_inv;

		// Check whether the program should be stopped early due to the li penalty
		
        pen(s + 1) = pen(s) - std::abs(gamma * A_A);

		if (!(std::isnan(l1)) && (((pen(s + 1) * 2) / l2_sqrt_inv) <= l1)) {
			
            double ps1_l2 = pen(s + 1) * 2 / l2_sqrt_inv;
			double 	ps_l2 = pen(s) * 2 / l2_sqrt_inv;
			beta_hat = (((ps_l2 - l1) / (ps_l2 - ps1_l2)) * beta_hat + ((l1 - ps1_l2) / (ps_l2 - ps1_l2)) * beta_curr) * l2_sqrt_inv;
			return beta_hat;
		
        }

		if (gamma_tilde < gamma_hat)
		{

			// Drop situation (for computational reasons cast the coefficients as zeros explicitly)
			
            drop = 1;
			beta_hat(gamma_index) = 0.;
			w_A_vec(gamma_index) = 0.;
			u_A.tail(N)(gamma_index) = 0.;
			++size_A_C;
			A_C_set(size_A_C - 1) = gamma_index;
			(A_set.array() == gamma_index).maxCoeff(&drop_index);

			if (drop_index != N - 1)
			{

				A_set(seq(drop_index, N - 2)) = A_set(seq(drop_index + 1, N - 1));
				sign(seq(drop_index, N - 2)) = sign(seq(drop_index + 1, N - 1));
			
            }

			A_set(N - 1) = -1;
			sign(N - 1) = 0;
			--size_A;

			// Downdate the Gram matrix
			
            CholUpdate(L, VectorXd::Zero(1), X(all, A_set.head(size_A)), drop_index, l2, "down");

			if (log)
			{
				if (curr_var_index != -1)
				{
					std::cout << '\n' << "Lars-step " << s << ": variable " << curr_var_index << " is added and " << gamma_index << " is dropped." << '\n';
				}
				else
				{
					std::cout << '\n' << "Lars-step " << s << ": variable " << gamma_index << " is dropped." << '\n';
				}
			}

		}
		else
		{

			// No drops

			if (N == size_A)
			{
				if (log) std::cout << '\n' << "Lars-step " << s << ": variable " << curr_var_index << " is added." << '\n';
				return beta_hat;
			}

			if (log)
			{
				if (curr_var_index != -1)
				{
					std::cout << '\n' << "Lars-step " << s << ": variable " << curr_var_index << " is added." << '\n';
				}
				else
				{
					std::cout << '\n' << "No variable has been added since previously a variable has been dropped." << '\n';
				}
			}
		}

		++s;
	}

	return beta_hat;
}

MatrixXd ElNetSolve::CD(
    const VectorXd& y, // VOI
    const MatrixXd& X, // Predictor Matriy
    const double& alpha, // KL1 and L2 mixing parameter
    double l1_start, // Starting value of the l1 grid
    double l1_end, // End value of the l1 grid
    double ratio, // ratio between the end and start value
    const int& grid_size, // Size of the l1 grind
    int selected, // Number of selected variables
    const double& comp_null, // Computational zero value
    const bool& log, // Talk to me (currently disabled)
    const int& max_iterations, // Maximum number of iteration steps for the inner and outer loop
    const double& conv_threshold // Conversion threshold
)
{

    // CD Algorithmus

    // Misshandling

    if (alpha < 0 || l1_start < 0 || grid_size < 0)
    {
        std::cout << '\n' << "Error! Either alpha, the starting lasso penalty, or the grid size is negative." << '\n';
        return VectorXd();
    }
    if (1 < alpha)
    {
        std::cout << '\n' << "Error! Alpha is bnound between 0 and 1." << '\n';
        return VectorXd();
    }

    // Dummies

    // Scalars

    int const N = X.cols(), T = X.rows();

    int active_set_size = 0;

    double l_max = (X.transpose() * y).cwiseAbs().maxCoeff() / (T * std::max(alpha, 0.001));
    l1_start = std::isnan(l1_start) ? l_max : l1_start;
    l1_end = std::isnan(l1_end) ? l1_start * std::max(ratio, 0.0001) : l1_end;
   
    // Create a sequence of lambdas

    if (l1_start <= l1_end)
    {
        std::cout << '\n' << "Error! l1_start = " << l1_start << " must be bigger then l_end = " << l1_end << "." << '\n';
        return EXIT_FAILURE * VectorXd::Ones(1);
    }

    VectorXd l1_seq = LogSeq(l1_start, l1_end, grid_size);
    MatrixXd beta_hat = MatrixXd::Zero(N, grid_size);

    // Pre-compute variables for computational purposes

    MatrixXd XtX = X.transpose() * X;
    VectorXd Xty = X.transpose() * y;
    double y_norm = y.squaredNorm();

    // loop over the l1 sequence

    for (int l = 0; l < grid_size; ++l) {

        // Use previous column as warm start

        VectorXd beta = beta_hat.col(std::max(0, l - 1));
        VectorXd beta_old = beta;
        bool converged = false;
        int iteration = 0;

        // Outer loop

        while (!converged && iteration < max_iterations) {
            
            converged = true;
            
            // Thresholding loop

            for (int n = 0; n < N; ++n) {
                
                // Compute estimate for given variable n

                double penalty = l1_seq(l) * alpha;
                double z = XtX(n, n) * beta(n) - (XtX.col(n).dot(beta) - XtX(n, n) * beta(n)) + Xty(n);
                double beta_n_old = beta(n);
                
                // Soft-thresholding

                if (z > penalty) {

                    beta(n) = (z - penalty) / (XtX(n, n) + l1_seq(l) * (1 - alpha));
                
                }
                else if (z < -penalty) {
                
                    beta(n) = (z + penalty) / (XtX(n, n) + l1_seq(l) * (1 - alpha));
                
                }
                else {
                
                    beta(n) = 0.0;
                
                }

                if (std::abs(beta(n) - beta_n_old) > conv_threshold) {

                    converged = false;

                }
            }

            iteration++;

            // Inner loop for active set estimate refinement
            
            if (!converged) {

                // Initialise the active set, i.e., the set of indices corresponding to non-zero estimates

                vector<int> active_set;

                for (int i = 0; i < N; ++i) {
                    if (beta(i) != 0) active_set.push_back(i);
                }

                active_set_size = active_set.size();

                bool active_set_converged = false;

                while (!active_set_converged && iteration < max_iterations) {

                    active_set_converged = true;

                    for (int i : active_set) {

                        // Compute the estimator updates

                        double penalty = l1_seq(l) * alpha;
                        double z = XtX(i, i) * beta(i) - (XtX.col(i).dot(beta) - XtX(i, i) * beta(i)) + Xty(i);
                        double beta_i_old = beta(i);

                        // Soft-thresholding

                        if (z > penalty) {

                            beta(i) = (z - penalty) / (XtX(i, i) + l1_seq(l) * (1 - alpha));
                        
                        }
                        else if (z < -penalty) {

                            beta(i) = (z + penalty) / (XtX(i, i) + l1_seq(l) * (1 - alpha));
                        
                        }
                        else {

                            beta(i) = 0.0;
                        
                        }

                        if (std::abs(beta(i) - beta_i_old) > conv_threshold) {
                            active_set_converged = false;
                        }

                    }

                    iteration++;
                }

            }

        }

        beta_hat.col(l) = beta;

        if (selected <= active_set_size) {

            return beta_hat;

        }

    }

    return beta_hat;

}