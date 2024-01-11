################################################################################
### GP Inference
################################################################################

#######
## Packages
#######

# Package names
packages <- c("fields","ggplot2", "dplyr","RandomFields")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

library(gpboost)

# Function for Simulating Data
source("https://raw.github.com/TimGyger/iterativeFSA/master/Data/Simulate_Data.R")

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

num_samples <- 10
l_mat <- list()
ind <- 1
for (ii in 1:num_samples) {
  set.seed(1)
  simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = ii)
  
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
  
  set.seed(100)
  simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = ii)
  
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
  coords_test <- as.matrix(X)
  y_test <- Y
  
  #######################
  ### Comparison
  #######################
  
  vec_info <- rep(0,8)
  # Initialize matrix
  for (i in 1:2) {
    if (i == 1){
      mim <- "cholesky"
    } else {
      mim <- "iterative"
    }
    t1 <- Sys.time()
    gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = mim,
                           likelihood = likelihood,seed = 10, num_ind_points = 500,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                           y = y_train,params = list(maxit=100,trace=TRUE,cg_delta_conv = 0.001,
                                                                             cg_preconditioner_type = "predictive_process_plus_diagonal",
                                                                             cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = 50,
                                                                             seed_rand_vec_trace = 10,reuse_rand_vec_trace = T,lr_cov = 1e-3))
    
    vec_info[4] <- Sys.time() - t1
    vec_info[1:3] <- gp_model$get_cov_pars()
    t1 <- Sys.time()
    gp_model$set_prediction_data(cg_delta_conv_pred = 1e-3, nsim_var_pred = 500)
    pred <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
    vec_info[5] <- sqrt(mean((y_test-pred$mu)^2))
    vec_info[6] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
    help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
    help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
    vec_info[7] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
    
    vec_info[8] <- Sys.time() - t1
    l_mat[[ind]] <- vec_info
    ind <- ind+1
  }
  
  
}

#######################
### Table
#######################
tab_final_chol <- l_mat[[1]]
tab_final_iter <- l_mat[[2]]
for (i in 2:num_samples) {
  tab_final_chol <- cbind(tab_final_chol,l_mat[[i*2-1]])
  tab_final_iter <- cbind(tab_final_iter,l_mat[[i*2]])
}
tab_final <- cbind(apply(tab_final_chol,1,mean),apply(tab_final_iter,1,mean))
tab_final_se <- cbind(apply(tab_final_chol,1,sd)/sqrt(num_samples),apply(tab_final_iter,1,sd)/sqrt(num_samples))

tab_final <- cbind(tab_final,tab_final_se)
colnames(tab_final) <- c("Exact","Iterative","Exact Standard Error","Iterative Standard Error")
rownames(tab_final) <- c("sigma","sigma_1","rho","time","RMSE","Log-Score","CRPSE","time_pred")


