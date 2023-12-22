################################################################################
### Comparison of Preconditioner (FITC-P vs. Pivoted Cholesky)
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
# number of data points
n <- 100000
likelihood <- "gaussian"


#######################
### Compare Runtime & CG-Iterations
#######################

# Effective ranges
vec_ER <- c(0.2,0.05,0.01)
# Initialize matrices and vectors
mat_time <- matrix(0,3,5)
colnames(mat_time) <- c("No Preconditioner","FITC-P","PIVCHOL 200","PIVCHOL 500","PIVCHOL 1000")
rownames(mat_time) <- c("EffRange 0.2","EffRange 0.01","EffRange 0.05")

for (i in vec_ER) {
  ii <- 1
  # Effective Range
  effectiverange <- i
  
  # Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
  arange <- effectiverange/4.7439
  
  ## Generate data points
  # Simulate
  set.seed(1)
  simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = 1)
  
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
  
  # Theoretically correct parameters
  if (likelihood == "gaussian") {
    init_cov_pars <- c(sigma_error,sigma,arange)
  } else {
    init_cov_pars <- c(sigma_error,arange)
  }
  
  for (j in 1:5) {
    
    if(j == 1){
      # Without Preconditioner
      num_neighbors_pred = 101
    } else if(j == 2){
      # FITC-P
      num_neighbors_pred = 100
    } else if(j == 3){
      # PivChol k = 200
      num_neighbors_pred = 200-1
    } else if(j == 4){
      # PivChol k = 500
      num_neighbors_pred = 500-1
    } else if(j == 5){
      # PivChol k = 1000
      num_neighbors_pred = 1000-1
    }
    
    t1 <- Sys.time()
    gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                        likelihood = likelihood,num_ind_points = 499,cov_fct_taper_range = 0.016, num_neighbors = 50,                       
                        gp_approx = "full_scale_tapering",num_neighbors_pred = num_neighbors_pred,
                        matrix_inversion_method = "cg",seed = 10)
        
    NEGLL <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
    mat_time[ii,j] <- Sys.time() - t1
    ii <- ii + 1
  }
}
