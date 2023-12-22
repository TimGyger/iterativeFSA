################################################################################
### GP Inference
################################################################################

#######
## Packages
#######

# Package names
packages <- c("fields","plotly", "dplyr","RandomFields","rlist","ggpubr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

library(gpboost)

# Function for Simulating Data
source("C:/Users/JumpStart/Desktop/Paper 1/Data/Simulate_Data.R")

#####################################################
# Parameters
#####################################################

## Covariance Function
# sigma
sigma <-  1
# Covariance Function
Covfunct <-  "Matern"
# Smoothness Parameter
smoothness <-  3/2
# Variance of error term
sigma_error = 1

## Data
n <- 100000
likelihood <- "gaussian"

# Effective Range
effectiverange <- 0.05

# Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
arange <- effectiverange/4.7439

set.seed(1)
simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = vec_ER[i],sigma_param = sigma, sigma_error = sigma_error,seed = 1)

# Y values
Y <- simdata[[1]]
# X values (locations)
X <- simdata[[2]]
# Actual number of data points
n <- dim(X)[1]

## Ordering
Y <- Y[order((X[,1])^2+(X[,2])^2)]
X <- X[order((X[,1])^2+(X[,2])^2),]

## Plot
quilt.plot(X[,1],X[,2],Y,nx = 200)
grid()
coords_train <- X
y_train <- Y


set.seed(10)
simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = vec_ER[i],sigma_param = sigma, sigma_error = sigma_error,seed = 1)

# Y values
Y <- simdata[[1]]
# X values (locations)
X <- simdata[[2]]
# Actual number of data points
n <- dim(X)[1]

## Ordering
Y <- Y[order((X[,1])^2+(X[,2])^2)]
X <- X[order((X[,1])^2+(X[,2])^2),]

## Plot
quilt.plot(X[,1],X[,2],Y,nx = 200)
grid()
coords_test <- X
y_test <- Y

#######################
### Comparison
#######################

l_mat <- list()
vec_info <- rep(0,8)
# Initialize matrix
for (i in 1:2) {
  if (i == 1){
    mim <- "cholesky"
  } else {
    mim <- "cg"
  }
  t1 <- Sys.time()
  gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = mim,
                         likelihood = likelihood,vecchia_ordering = "random",seed = 10, num_ind_points = 499,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                         y = y_train,num_neighbors_pred = 6, params = list(maxit=1000,trace=TRUE,cg_delta_conv = 0.001,
                                                                            cg_preconditioner_type = "predictive_process_plus_diagonal", piv_chol_rank = 500,
                                                                            cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = 50,
                                                                            seed_rand_vec_trace = 10,reuse_rand_vec_trace = T,lr_cov = 1e-8))
  
  vec_info[4] <- Sys.time() - t1
  vec_info[1:3] <- gp_model$get_cov_pars()
  t1 <- Sys.time()
  pred <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
  vec_info[8] <- Sys.time() - t1
  l_mat <- list.append(l_mat,vec_info)
}
