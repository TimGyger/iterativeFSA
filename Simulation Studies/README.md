# Simulation Studies

Each script conducts a specific simulation study and includes a plotting function at the end to visualize the results.
The simulated data is generated using the ```sim_data()``` function, located in ```root\Data\Simulation\Simulate_Data.R```.

## Methods for choosing inducing points

The script ```Inducing_Point_Method.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and different inducing point methods (random selection, covertree, kmeans++) for the FITC and full-scale approximation and generates Figure 1 in the main manuscript in Section ```Methods for choosing inducing points``` and Figure 1 in the Appendix ```Additional plots and tables```. 
On a standard laptop, memory limitations may arise; therefore, it is recommended to use no more than 500 inducing points to prevent such issues.

These computations can be time-intensive on typical hardware, especially for the full-scale approximation. Using Cholesky-based methods with the setting ```matrix_inversion_method = "iterative"``` can significantly accelerate the process in this case. 
However, note that this approach may yield a less accurate estimation of the negative log-likelihood.

## Choice of full-scale approximation parameters

The script ```Model_Parameter.R``` calculates the negative log-likelihood at the data-generating parameters for various numbers of inducing points and taper ranges and generates Figure 1 in the Appendix ```Choice of full-scale approximation parameters```. 
On a standard laptop, memory limitations may arise; therefore, it is recommended to use no more than 500 inducing points to prevent such issues.

These computations can be time-intensive on typical hardware. Using Cholesky-based methods with the setting ```matrix_inversion_method = "iterative"``` can significantly accelerate the process. 
However, note that this approach may yield a less accurate estimation of the negative log-likelihood.
