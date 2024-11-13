################################################################################
### Variance Analysis (NEGLL)
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
n <- 8000
likelihood <- "gaussian"


#######################
### NEGLL
#######################


# Effective ranges
vec_ER <- c(0.5,0.2,0.05)
# Number of sample vectors
vec_samples <- c(5,10,30,50,80,100)
# Initialize list
l_mat <- list()
ind <- 1
for (i in 1:3) {
  i <- 1
  # Effective Range
  effectiverange <- vec_ER[i]
  
  # Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
  arange <- effectiverange/4.7439
  
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
  a <- Sys.time()
  gp_model1 <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                         likelihood = likelihood,seed = 10, num_ind_points = 50,cov_fct_taper_range = 2,gp_approx = "vecchia",
                         vecchia_ordering = "none",num_neighbors =20,
                         y = y_train,params = list(maxit=8,trace=TRUE,optimizer_cov = "gradient_descent",
                                                   init_cov_pars = c(0.5,0.7,0.5)))
  
  set.seed(2)
  simdata <- sim_data(n = 10000,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error,seed = 1)
  
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
  set_prediction_data(gp_model1,vecchia_pred_type = "order_obs_first_cond_obs_only")
  pred21 <- predict(gp_model1,gp_coords_pred = coords_test, y = as.numeric(y_train), predict_var = T)
  Sys.time()-a
  NEGLL <- gp_model$get_current_neg_log_likelihood()
  mat1 <- matrix(0,length(vec_samples),25)
  mat2 <- matrix(0,length(vec_samples),25)
  for (ii in 1:length(vec_samples)) {
    for (j in 1:2){
      vec_s <- rep(0,25)
      for (jj in 1:25){
        
        if (j == 2){
          gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "iterative",
                                 likelihood = likelihood,seed = 10, num_ind_points = 500,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                                 y = y_train,params = list(maxit=0,trace=TRUE,cg_delta_conv = 0.001,optimizer_cov = "gradient_descent",
                                                           cg_preconditioner_type = "predictive_process_plus_diagonal",
                                                           cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = vec_samples[ii],
                                                           seed_rand_vec_trace = jj*10,reuse_rand_vec_trace = T,init_cov_pars = c(1,1,arange)))
          mat1[ii,jj] <- (gp_model$get_current_neg_log_likelihood()-NEGLL)/NEGLL
        } else {
          gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "iterative",
                                 likelihood = likelihood,seed = 10, num_ind_points = 500,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                                 y = y_train,params = list(maxit=0,trace=TRUE,cg_delta_conv = 0.001,optimizer_cov = "gradient_descent",
                                                           cg_preconditioner_type = "none",
                                                           cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = vec_samples[ii],
                                                           seed_rand_vec_trace = jj*10,reuse_rand_vec_trace = T,init_cov_pars = c(1,1,arange)))
          mat2[ii,jj] <- (gp_model$get_current_neg_log_likelihood()-NEGLL)/NEGLL
        }
        
      }
      print(mat1)
      print(mat2)
    }
  }
  l_mat[[ind]] <- mat1
  l_mat[[ind+1]] <- mat2
  ind <- ind+2
}


###################
### Plots
###################

l_plots <- vector('list', 6)
for (i in 1:6) {
  mat_NEGLL <- l_mat[[i]]
  mat_NEGLL2 <- cbind(vec_samples,mat_NEGLL)
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat)[1] <- c("x")
  data_long <- reshape2::melt(data_mat, id = "x")
  
  l_plots[[i]] <- local({
    ggplot(data_long,            
           aes(x = x,
               y = value,
               color = variable)) +  geom_line() + geom_point() +  scale_colour_manual(values = c(2:3),name = "", 
                                                                                       labels = c("FITC-P","No Preconditioner"))+
      theme(
        legend.position="bottom",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = "Number of Sample Vectors", y = "") +
      scale_x_continuous(breaks =vec_samples) + scale_y_log10()
  })

}
ggarrange(l_plots[[1]], l_plots[[2]], l_plots[[3]],
          l_plots[[4]], l_plots[[5]], l_plots[[6]],ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
