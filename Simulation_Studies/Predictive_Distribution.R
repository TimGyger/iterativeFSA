################################################################################
### Comparison of Predictive Distribution (Not possible to replicate for Lanczos)
################################################################################

#######
## Packages
#######

source("https://raw.githubusercontent.com/TimGyger/iterativeFSA/refs/heads/main/Packages.R")

#######
## Data
#######

source("https://raw.github.com/TimGyger/iterativeFSA/master/Data/Simulation/Simulate_Data.R")

#####################################################
# Parameters
#####################################################

### Toy example
Toy <- TRUE

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
if(Toy){
  n <- 10000
}
likelihood <- "gaussian"

# Effective Range
effectiverange <- 0.2

# Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
arange <- effectiverange/2.7

set.seed(1)
simdata <- sim_data(n = 2*n,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = ii)
  
# Y values
Y <- simdata[[1]][1:n]
# X values (locations)
X <- simdata[[2]][1:n,]
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
  
# Y values
Y <- simdata[[1]][(n+1):(2*n)]
# X values (locations)
X <- simdata[[2]][(n+1):(2*n),]
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
### Comparison Prediction
#######################


# Number of sample vectors / Rank
vec_samples <- c(50,200,500,1000,2000,5000)
if(Toy){
  vec_samples <- c(5,10,20,30,40,50)
}
# Initialize matrix
mat <- matrix(0,24,2)
rownames(mat) <- c("Log-Score 50","CRPS 50","RMSE 50","Time 50",
                   "Log-Score 200","CRPS 200","RMSE 200","Time 200",
                   "Log-Score 500","CRPS 500","RMSE 500","Time 500",
                   "Log-Score 1000","CRPS 1000","RMSE 1000","Time 1000",
                   "Log-Score 2000","CRPS 2000","RMSE 2000","Time 2000",
                   "Log-Score 5000","CRPS 5000","RMSE 5000","Time 5000")
if(Toy){
  rownames(mat) <- c("Log-Score 5","CRPS 5","RMSE 5","Time 5",
                     "Log-Score 10","CRPS 10","RMSE 10","Time 10",
                     "Log-Score 20","CRPS 20","RMSE 20","Time 20",
                     "Log-Score 30","CRPS 30","RMSE 30","Time 30",
                     "Log-Score 40","CRPS 40","RMSE 40","Time 40",
                     "Log-Score 50","CRPS 50","RMSE 50","Time 50")
}
colnames(mat) <- c("Stochastic","Cholesky")
mat_var <- mat
gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = likelihood,seed = 10, num_ind_points = 500,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                       y = y_train,params = list(maxit=0,trace=TRUE,lr_cov = 1e-8,init_cov_pars = c(1,1,arange)))
t1 <- proc.time()[[3]]
pred_Chol <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
mat[4,2] <- proc.time()[[3]] - t1 
mat[1,2] <- -mean(dnorm(y_test,pred_Chol$mu,sqrt(pred_Chol$var), log = T))
help_vec <- dnorm((y_test-pred_Chol$mu)/sqrt(pred_Chol$var),rep(0,length(pred_Chol$mu)),rep(1,length(pred_Chol$mu)))
help_vec2 <- pnorm((y_test-pred_Chol$mu)/sqrt(pred_Chol$var),rep(0,length(pred_Chol$mu)),rep(1,length(pred_Chol$mu)))
mat[2,2] <- -mean(sqrt(pred_Chol$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred_Chol$mu)/sqrt(pred_Chol$var))*(2*help_vec2-1)))
k <- 0
for (ii in 1:length(vec_samples)) {
  num_it <- 25
  if(Toy){
    num_it <- 5
  }
  mat_rep <- matrix(0,4,num_it)
  for (jj in 1:num_it) {
    num_ind_points <- 500
    cov_fct_taper_range <- 0.016
    if(Toy){
      num_ind_points <- 200
      cov_fct_taper_range <- 0.052
    }
    gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "iterative",
                           likelihood = likelihood,seed = 10*jj, num_ind_points = num_ind_points,cov_fct_taper_range = cov_fct_taper_range,gp_approx = "full_scale_tapering",
                           y = y_train,params = list(maxit=0,trace=TRUE,cg_delta_conv = 0.001,init_cov_pars = c(1,1,arange),
                                                     cg_preconditioner_type = "predictive_process_plus_diagonal",
                                                     cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = 50,
                                                     seed_rand_vec_trace = 10,reuse_rand_vec_trace = T,lr_cov = 1e-8))
    t1 <- proc.time()[[3]]
    gp_model$set_prediction_data(nsim_var_pred = vec_samples[ii])
    pred <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
    mat_rep[1,jj] <- proc.time()[[3]] - t1
    mat_rep[2,jj] <- sqrt(mean((pred$var-as.numeric(unlist(pred_Chol$var)))^2))
    mat_rep[3,jj] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
    help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
    help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
    mat_rep[4,jj] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
  }
  mat[k+4,1] <- mean(mat_rep[1,])
  mat[k+3,1] <- mean(mat_rep[2,])
  mat[k+1,1] <- mean(mat_rep[3,])
  mat[k+2,1] <- mean(mat_rep[4,])
  mat_var[k+4,1] <- var(mat_rep[1,])
  mat_var[k+3,1] <- var(mat_rep[2,])
  mat_var[k+1,1] <- var(mat_rep[3,])
  mat_var[k+2,1] <- var(mat_rep[4,])
  k <- k + 4
}

print(mat)
print(mat_var)
