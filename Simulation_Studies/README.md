# Simulation Studies

Each script conducts a specific simulation study and includes a plotting function at the end to visualize the results.
The simulated data is generated using the ```sim_data()``` function, located in ```root\Data\Simulation\Simulate_Data.R```.

## Methods for choosing inducing points

The script ```Inducing_Point_Method.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and different inducing point methods (random selection, covertree, kmeans++) for the FITC and full-scale approximation and generates Figure 1 in the main manuscript in Section ```Methods for choosing inducing points``` and Figure 1 in the Appendix ```Additional plots and tables```. 
On a standard laptop, memory limitations may arise; therefore, it is recommended to use no more than 500 inducing points to prevent such issues.

These computations can be time-intensive on typical hardware, especially for the full-scale approximation. Using Cholesky-based methods with the setting ```matrix_inversion_method = "iterative"``` can significantly accelerate the process in this case. 
However, note that this approach may yield a less accurate estimation of the negative log-likelihood.

## Runtime Comparison of iterative and Cholesky-based methods

The script ```Comparison_Iterative_vs_Cholesky.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points, taper ranges and sample sizes and generates Figure 5 in the main manuscript in Section ```Accuracy and computational time of parameter estimates and
predictive distributions```. 
On a standard laptop, memory limitations may arise; therefore, it is recommended to use no more than 500 inducing points and n <= 100'000 to prevent such issues.

## Choice of full-scale approximation parameters

The script ```Model_Parameter.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and taper ranges and generates Figure SM1 in the Appendix ```Choice of full-scale approximation parameters```. 
On a standard laptop, memory limitations may arise; therefore, it is recommended to use no more than 500 inducing points to prevent such issues.

These computations can be time-intensive on typical hardware. Using Cholesky-based methods with the setting ```matrix_inversion_method = "iterative"``` can significantly accelerate the process. 
However, note that this approach may yield a less accurate estimation of the negative log-likelihood.

## Variance Analysis

The script ```Variance_Analysis_NEGLL.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of sample vectors and generates box-plots. For the negative log-likelihood at the initial parameters remove ```initial_cov_pars = ...```. It generates Figures SM3 - SM5 in the Appendix.

The script ```Variance_Analysis_Grad_NEGLL.R``` calculates the gradients of the negative log-likelihood at the data-generating parameters for various numbers of sample vectors and generates box-plots. It generates Figures SM6 - SM8 in the Appendix.

Set ```n = 1000``` for faster computations.

## Predictive Distribution

The script ```Predictive_Distribution.R``` calculates the predictive mean and variances for different number of sample vectors and generates the table used for
Figure 3 in the main mansucript.

Set ```n = 1000``` for faster computations.

## Estimation and Prediction

The script ```GP_Inference.R``` estimates the parameters and calculates the predictive mean and variances for different number of sample vectors and generates the table used for Figure 4 in the main mansucript.

Set ```n = 1000``` for faster computations.

## Preconditioner Comparison for Vecchia approximation

The script ```Comparison_Preconditioner_Vecchia.R``` compares FITC, pivoted Cholesky, and Vecchia approximation with diagonal update (VADU) preconditioners in terms of runtime and the variance of marginal likelihood estimates. Set ```n = 1000``` for faster computations. It generates Figure 7 in the main manuscript.

The script ```FITC_vs_pivotedCholesky_Vecchia.R``` compares FITC and pivoted Cholesky preconditioners in terms of runtime and the variance of marginal likelihood estimates for various ranks and number of inducing points. Set ```n = 1000``` for faster computations. It generates Figure 6 in the main manuscript.


