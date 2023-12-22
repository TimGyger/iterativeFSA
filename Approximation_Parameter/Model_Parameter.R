################################################################################
### Analysis of approximation parameter
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
### Compare exact NEGLL
#######################

# Effective ranges
vec_ER <- c(0.2,0.05,0.01)
# Number of inducing points
vec_ind_points <- c(100,200,300,400,500,600,700,800,900,1000)-1
# Taper Range
vec_taper_range <- c(0.0055,0.0079,0.0098,0.0112,0.0126,0.0138,0.015,0.016,0.017,0.0179,0.0188,0.0197)

# Initialize matrices and vectors
mat_NEGLL <- matrix(0,12,10)
colnames(mat_NEGLL) <- vec_ind_points
l_ap <- list()


for (i in vec_ER) {
  
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
  j_ind <- 1
  for (j in vec_ind_points) {
    i_ind <- 1
    for (ji in vec_taper_range) {
     
      gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,
                          likelihood = likelihood,num_ind_points = j,                     
                          gp_approx = "full_scale_tapering", cov_fct_shape = 1.5,cov_fct_taper_range = ji,  
                          matrix_inversion_method = "cholesky",seed = 10)
      
      mat_NEGLL[i_ind,j_ind] <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
      i_ind <- i_ind + 1
    }
    j_ind <- j_ind + 1
    
  }
  l_ap <- list.append(mat_NEGLL)
}


###################
### Plots
###################

FITC_plots <- vector('list', 3)
for (i in 1:3) {
  
  mat_NEGLL1 <- l_ap[[i]]
  mat_NEGLL2 <- cbind(c(1:12)*10,mat_NEGLL)
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat)[1] <- c("x")
  data_long <- reshape2::melt(data_mat, id = "x")
  
  AP_plots[[i]] <- local({
    ggplot(data_long,            
           aes(x = x,
               y = value,
               color = variable)) +  geom_line() + geom_point() +  scale_colour_grey(start = 0, end = 0.8,name = "Number of inducing points", 
                                                                                     labels = c("100","200","300","400","500","600","700","800",
                                                                                                "900","1000"))+
      theme(
        legend.position="bottom",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = "Average number of non-zeros per row", y = "NEGLL") +
      scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100,110,120))
  })
  
  ggarrange(AP_plots[[1]], AP_plots[[2]], AP_plots[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
}




