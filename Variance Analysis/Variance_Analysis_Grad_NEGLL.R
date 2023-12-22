################################################################################
### Variance Analysis (NEGLL)
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

#######################
### Gradient of NEGLL
#######################


# Effective ranges
vec_ER <- c(0.2,0.05,0.01)
# Number of sample vectors
vec_samples <- c(5,10,30,50,80,100)
# Initialize list
l_mat <- list()

Var_Mat_relerr <- matrix(0,length(vec_samples),3)
Range_Mat_relerr <- matrix(0,length(vec_samples),3)
Var_Mat_var <- matrix(0,length(vec_samples),3)
Range_Mat_var <- matrix(0,length(vec_samples),3)

gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = likelihood,seed = 10, num_ind_points = 499,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                       y = y_train,params = list(maxit=1,trace=TRUE,lr_cov = 1e-8))

Var_para <- (log(1)-log(gp_model$get_cov_pars()[2]/gp_model$get_cov_pars()[1]))/1e-8
Range_para <- (log(gp_model$get_optim_params()$init_cov_pars[3])-log(gp_model$get_cov_pars()[3]))/-1e-8

for (i in 1:3) {
  prec <- "predictive_process_plus_diagonal"
  cc <- "none"
  if(i == 2){
    cc <- "random"
  }
  if(i == 3){
    prec <- "none"
  }
  for (jj in 1:length(vec_num_samples)) {
    vec_s <- rep(0,num_trials)
    vec_s1 <- rep(0,num_trials)
    for (ii in 1:num_trials) {
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cg",
                             likelihood = likelihood,vecchia_ordering = cc,seed = 10, num_ind_points = 499,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                             y = y_train,params = list(maxit=1,trace=TRUE,cg_delta_conv = 0.001,
                                                       cg_preconditioner_type = prec,
                                                       cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = jj,
                                                       seed_rand_vec_trace = ii*10,reuse_rand_vec_trace = T,lr_cov = 1e-8))
      
      vec_s[ii] <- abs((log(1)-log(gp_model$get_cov_pars()[2]/gp_model$get_cov_pars()[1]))/1e-8-Var_para)/Var_para
      vec_s1[ii] <- abs((log(gp_model$get_optim_params()$init_cov_pars[3])-log(gp_model$get_cov_pars()[3]))/-1e-8-Range_para)/Range_para
    }
    Var_Mat_relerr[jj,i] <- mean(vec_s)
    Range_Mat_relerr[jj,i] <- mean(vec_s1)
    Var_Mat_var[jj,i] <- var(vec_s)
    Range_Mat_var[jj,i] <- var(vec_s1)
  }
}
l_mat <- list.append(l_mat,Var_Mat_relerr,Var_Mat_var,Range_Mat_relerr,Range_Mat_var)

###################
### Plots
###################

l_plots <- vector('list', 4)
for (i in 1:4) {
  mat_NEGLL1 <- l_mat[[i]]
  mat_NEGLL2 <- cbind(vec,mat_NEGLL)
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat)[1] <- c("x")
  data_long <- reshape2::melt(data_mat, id = "x")
  
  l_plots[[i]] <- local({
    ggplot(data_long,            
           aes(x = x,
               y = value,
               color = variable)) +  geom_line() + geom_point() +  scale_colour_manual(name = "Number of Sample Vectors", 
                                                                                       labels = c("5","10","50","80","100"))+
      theme(
        legend.position="bottom",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = "Number of Sample Vectors", y = "") +
      scale_x_continuous(breaks =vec1)
  })
  
  ggarrange(l_plots[[1]], l_plots[[2]], l_plots[[3]],
            l_plots[[4]], ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
}
