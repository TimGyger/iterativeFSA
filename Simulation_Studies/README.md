# Simulation Studies

Each script conducts a specific simulation study and includes a plotting function at the end to visualize the results.
The simulated data is generated using the ```sim_data()``` function, located in ```root\Data\Simulation\Simulate_Data.R```.

## Toy Example

Each script includes the variable ```Toy```, which activates a simplified "Toy example" mode when set to ```True```. This setting significantly reduces computation time and memory requirements by using a smaller sample size and scaling down other parameters. However, in this Toy mode, some plots may differ from those presented in the manuscript. For instance, runtimes for Cholesky-based computations are highly dependent on the sample size, and these differences become apparent in the Toy setting.

## Methods for choosing inducing points

The script ```Inducing_Point_Method.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and different inducing point methods (random selection, covertree, kmeans++) for the FITC and full-scale approximation and generates Figure 1 in the main manuscript in Section ```Methods for choosing inducing points``` and Figure 1 in the Appendix ```Additional plots and tables```. 

## Runtime Comparison of iterative and Cholesky-based methods

The script ```Comparison_Iterative_vs_Cholesky.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points, taper ranges and sample sizes and generates Figure 4 in the main manuscript in Section ```Computational time and accuracy of log-marginal likelihoods, parameter estimates, and predictive distributions```. 

## Choice of full-scale approximation parameters

The script ```Model_Parameter.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and taper ranges and generates Figure SM1 in the Appendix ```Choice of full-scale approximation parameters```. 

## Variance Analysis

The script ```Variance_Analysis_NEGLL.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of sample vectors and generates box-plots. For the negative log-likelihood at the initial parameters remove ```initial_cov_pars = ...```. It generates Figures SM3 - SM5 in the Appendix.

The script ```Variance_Analysis_Grad_NEGLL.R``` calculates the gradients of the negative log-likelihood at the data-generating parameters for various numbers of sample vectors and generates box-plots. It generates Figures SM6 - SM8 in the Appendix.

## Predictive Distribution

The script ```Predictive_Distribution.R``` calculates the predictive mean and variances for different number of sample vectors and generates the table used for
Figure 3 in the main mansucript.

## Estimation and Prediction

The script ```GP_Inference.R``` estimates the parameters and calculates the predictive mean and variances for different number of sample vectors and generates the table used for Figure 5 in the main mansucript.

## Preconditioner Comparison for Vecchia approximation

The script ```Comparison_Preconditioner_Vecchia.R``` compares FITC, pivoted Cholesky, and Vecchia approximation with diagonal update (VADU) preconditioners in terms of runtime and the variance of marginal likelihood estimates. It generates Figure 6 in the main manuscript.

The script ```FITC_vs_pivotedCholesky_Vecchia.R``` compares FITC and pivoted Cholesky preconditioners in terms of runtime and the variance of marginal likelihood estimates for various ranks and number of inducing points. It generates Figure SM9 in the Appendix.


