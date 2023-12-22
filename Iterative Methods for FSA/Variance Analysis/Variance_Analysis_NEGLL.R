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


#######################
### NEGLL
#######################


# Effective ranges
vec_ER <- c(0.2,0.05,0.01)
# Number of sample vectors
vec_samples <- c(5,10,30,50,80,100)
# Initialize list
l_mat <- list()

for (i in 1:3) {
  
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
  
  gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                      likelihood = likelihood,num_ind_points = 499,cov_fct_taper_range = 0.016, num_neighbors = 50,                       
                      gp_approx = "full_scale_tapering",
                      matrix_inversion_method = "cholesky",seed = 10)
  
  NEGLL <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
  mat1 <- matrix(0,length(vec_samples),2)
  mat2 <- matrix(0,length(vec_samples),2)
  for (ii in 1:length(vec_samples)) {
    for (j in 1:2){
      if (j == 1){
        num_neighbors_pred = 101
      } else {
        num_neighbors_pred = 100
      }
      vec_s <- rep(0,25)
      for (jj in 1:25){
        gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                            likelihood = likelihood,num_ind_points = 499,cov_fct_taper_range = 0.016, num_neighbors = vec_samples[ii],                       
                            gp_approx = "full_scale_tapering",num_neighbors_pred = num_neighbors_pred,
                            matrix_inversion_method = "cg",seed = jj*10)
        
        vec_s[jj] <- abs(gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)-NEGLL)/NEGLL
      }
      mat1[ii,j] <- mean(vec_s)
      mat2[ii,j] <- var(vec_s)
    }
  }
  l_mat <- list.append(l_mat,mat1,mat2)
}


###################
### Plots
###################

l_plots <- vector('list', 6)
for (i in 1:6) {
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
            l_plots[[4]], l_plots[[5]], l_plots[[6]],ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
}
