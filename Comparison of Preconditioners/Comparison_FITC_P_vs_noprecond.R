################################################################################
### Comparison of FITC-P vs. No-Preconditioner vs. exact
################################################################################

#######
## Packages
#######

# Package names
packages <- c("fields","ggplot2", "dplyr","RandomFields","ggpubr")

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
effectiverange <- 0.2

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

#######################
### Compare Runtime and CG-Iterations
#######################


# Effective ranges
vec_ER <- c(0.2,0.05,0.01)
# Number of inducing points
vec_ind_points <- c(100,200,300,400,500,600,700,800,900,1000)
# Taper Range
vec_taper_range <- c(0.0055,0.0079,0.0098,0.0112,0.0126,0.0138,0.015,0.016,0.017,0.0179,0.0188,0.0197)
# sample size
n_vec <- c(10,20,50,80,100,150,200)*1000
# Taper range 2
vec_taper_range2 <- c(0.05,0.036,0.0227,0.0179,0.016,0.013,0.0113)
# Initialize list
l_mat <- list()
ind <- 1

for (i in 1:3) {
  # Inducing Points
  if (i == 1){
    mat_time <- matrix(0,length(vec_ind_points) ,3)
    iii <- 1
    for (ii in vec_ind_points) {
      
      for (jj in 1:3) {
        if(jj == 1){
          mm <- "cholesky"
        } else {
          mm <- "iterative"
        }
        t1 <- Sys.time()
        gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                            likelihood = likelihood,num_ind_points = ii,cov_fct_taper_range = 0.016,                     
                            gp_approx = "full_scale_tapering",
                            matrix_inversion_method = mm,seed = 10)
        if (jj == 3){
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "predictive_process_plus_diagonal"))
        } else {
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "none"))
        }
        NEGLL <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
        mat_time[iii,jj] <- Sys.time() - t1
      }
      iii <- iii + 1
    }
    l_mat[[ind]] <- mat_time
    ind <- ind+1
  }
  
  if (i == 2){
    mat_time <- matrix(0,length(vec_taper_range) ,3)
    iii <- 1
    for (ii in vec_taper_range) {
      
      for (jj in 1:3) {
        if(jj == 1){
          mm <- "cholesky"
        } else {
          mm <- "iterative"
        }
        t1 <- Sys.time()
        gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                            likelihood = likelihood,num_ind_points = 500,cov_fct_taper_range = ii,                   
                            gp_approx = "full_scale_tapering",
                            matrix_inversion_method = mm,seed = 10)
        if (jj == 3){
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "predictive_process_plus_diagonal"))
        } else {
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "none"))
        }
        NEGLL <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
        mat_time[iii,jj] <- Sys.time() - t1
      }
      iii <- iii + 1
    }
    l_mat[[ind]] <- mat_time
    ind <- ind+1
  }
  
  if (i == 3){
    mat_time <- matrix(0,length(n_vec) ,3)
    iii <- 1
    for (ii in n_vec) {
      set.seed(1)
      simdata <- sim_data(n = ii,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = 1)
      
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
      
      for (jj in 1:3) {
        if(jj == 1){
          mm <- "cholesky"
        } else {
          mm <- "iterative"
        }
        t1 <- Sys.time()
        gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                            likelihood = likelihood,num_ind_points = 500,cov_fct_taper_range = vec_taper_range2[iii],                       
                            gp_approx = "full_scale_tapering",
                            matrix_inversion_method = mm,seed = 10)
        if (jj == 3){
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "predictive_process_plus_diagonal"))
        } else {
          gp_model$set_optim_params(params = list(cg_preconditioner_type = "none"))
        }
        NEGLL <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
        mat_time[iii,jj] <- Sys.time() - t1
      }
      iii <- iii + 1
    }
    l_mat[[ind]] <- mat_time
    ind <- ind+1
  }
}


###################
### Plots
###################

l_plots <- vector('list', 3)
for (i in 1:3) {
  if(i == 1){
    xxlab <- "Number of inducing points"
    vec <- vec_ind_points + 1
    vec1 <- vec
  } else if(i == 2){
    xxlab <- "Average number of nonzero entries per row"
    vec <- vec_taper_range
    vec1 <- vec
  } else {
    xxlab <- "Sample size in thousands"
    vec <- n_vec
    vec1 <- vec/1000
  }
  mat_NEGLL1 <- l_mat[[i]]
  mat_NEGLL2 <- cbind(vec,mat_NEGLL1)
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat)[1] <- c("x")
  data_long <- reshape2::melt(data_mat, id = "x")
  
  l_plots[[i]] <- local({
    ggplot(data_long,            
           aes(x = x,
               y = value,
               color = variable)) +  geom_line() + geom_point() +  scale_colour_manual(values = c(2:4),name = "", 
                                                                                       labels = c("Cholesky","No Preconditioner",
                                                                                                  "FITC-P"))+
      theme(
        legend.position="bottom",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = xxlab, y = "Runtime in seconds (s)") +
      scale_x_continuous(breaks =vec1)
  })
  
}

ggarrange(l_plots[[1]], l_plots[[2]], l_plots[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
