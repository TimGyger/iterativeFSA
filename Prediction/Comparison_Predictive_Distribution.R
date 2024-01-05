################################################################################
### Comparison of Predictive Distribution (Not possible to replicate for Lanczos)
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
effectiverange <- 0.05

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


set.seed(100)
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
coords_test <- X
y_test <- Y

#######################
### Comparison Prediction
#######################


# Number of sample vectors / Rank
vec_samples <- c(50,100,200,500,1000,5000)
# Initialize matrix
mat <- matrix(0,24,3)
rownames(mat) <- c("Log-Score 50","CRPS 50","RMSE 50","Time 50",
                   "Log-Score 100","CRPS 100","RMSE 100","Time 100",
                   "Log-Score 200","CRPS 200","RMSE 200","Time 200",
                   "Log-Score 500","CRPS 500","RMSE 500","Time 500",
                   "Log-Score 1000","CRPS 1000","RMSE 1000","Time 1000",
                   "Log-Score 5000","CRPS 5000","RMSE 5000","Time 5000")
colnames(mat) <- c("Lanczos","Stochastic","Exact")
mat_var <- mat
gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = likelihood,seed = 10, num_ind_points = 499,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                       y = y_train,params = list(maxit=0,trace=TRUE,lr_cov = 1e-8,init_cov_pars = init_cov_pars))
t1 <- Sys.time()
pred_Chol <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
mat[4,3] <- Sys.time() - t1 
mat[1,3] <- -mean(dnorm(y_test,pred_Chol$mu,sqrt(pred_Chol$var), log = T))
help_vec <- dnorm((y_test-pred_Chol$mu)/sqrt(pred_Chol$var),rep(0,length(pred_Chol$mu)),rep(1,length(pred_Chol$mu)))
help_vec2 <- pnorm((y_test-pred_Chol$mu)/sqrt(pred_Chol$var),rep(0,length(pred_Chol$mu)),rep(1,length(pred_Chol$mu)))
mat[2,3] <- -mean(sqrt(pred_Chol$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred_Chol$mu)/sqrt(pred_Chol$var))*(2*help_vec2-1)))
for (i in 1:2) {
  k <- 0
  if (i == 1){
    nn <- 2
  } else {
    nn <- 6
  }
  for (ii in 1:length(vec_samples)) {
    mat_rep <- matrix(0,4,25)
    for (jj in 1:25) {
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cg",
                             likelihood = likelihood,vecchia_ordering = "random",seed = 10, num_ind_points = 499,cov_fct_taper_range = 0.016,gp_approx = "full_scale_tapering",
                             y = y_train,num_neighbors_pred = nn, params = list(maxit=1,trace=TRUE,cg_delta_conv = 0.001,
                                                                                cg_preconditioner_type = "predictive_process_plus_diagonal", piv_chol_rank = vec_samples[ii],
                                                                                cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = 50,
                                                                                seed_rand_vec_trace = 10*jj,reuse_rand_vec_trace = T,lr_cov = 1e-8))
      t1 <- Sys.time()
      pred <- predict(gp_model, gp_coords_pred = coords_test, predict_var = T,y = y_train)
      mat_rep[1,jj] <- Sys.time() - t1
      mat_rep[2,jj] <- sqrt(mean((pred$var-pred_Chol$var)^2))
      mat_rep[3,jj] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
      help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
      help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
      mat_rep[4,jj] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
    }
    
    mat[k+4,i] <- mean(mat_rep[1,])
    mat[k+3,i] <- mean(mat_rep[2,])
    mat[k+1,i] <- mean(mat_rep[3,])
    mat[k+2,i] <- mean(mat_rep[4,])
    mat_var[k+4,i] <- var(mat_rep[1,])
    mat_var[k+3,i] <- var(mat_rep[2,])
    mat_var[k+1,i] <- var(mat_rep[3,])
    mat_var[k+2,i] <- var(mat_rep[4,])
    k <- k + 4
  }
}


###################
### Variance Table
###################

print(mat_var)

###################
### Plots
###################

vec_num1 <- vec_samples
vec_num1[-2] <- ""
vec_num2 <- vec_samples
vec_num2[2] <- ""

df3 <- cbind(mat[c(1:6)*4-1,2],mat[c(1:6)*4-3,2],mat[c(1:6)*4,2],mat[c(1:6)*4-1,1],mat[c(1:6)*4-3,1],mat[c(1:6)*4,1])
colnames(df3) <- c("RMSE1","LogScore1","Time1","RMSE","LogScore","Time")
df3 <- as.data.frame(df3)
gfg_plot1 <- ggplot(data = df3) +  geom_line(aes(x = Time1,y = LogScore1,color = "Lanczos")) + 
  geom_point(aes(x = Time1,y = LogScore1), color = 2)  + 
  geom_text(aes(x = Time1,y = LogScore1, label = vec_num2),hjust=0.5, vjust=-0.5,color = 2,size = 4) +
  geom_text(aes(x = Time1,y = LogScore1, label = vec_num1),hjust=0.5, vjust=1.5,color = 2,size = 4) +
  geom_line(aes(x = Time,y = LogScore, color = "Stochastic")) + geom_point(aes(x = Time,y = LogScore), color = 4) +
  geom_text(aes(x = Time,y = LogScore, label = vec_num2),hjust=0.5, vjust=-0.5,color = 4,size = 4) +
  geom_text(aes(x = Time,y = LogScore, label = vec_num1),hjust=0.5, vjust=1.5,color = 4,size = 4) +
  geom_text(aes(x = 1400,y = mat[1,3], label = paste0("Exact: ",mat[4,3]," s")),size = 5) +
  scale_x_log10(limits = c(55,2000), breaks = c(50,100,200,500,1000,2000)) + scale_y_log10() + 
  labs(x = "Time (s)", y = "Log-Score") + theme(
    legend.position="bottom",
    legend.title=element_blank(),
    legend.background = element_rect(linetype="solid", 
                                     colour ="black")) + theme(text = element_text(size=20)) + 
  geom_hline(yintercept=mat[1,3], linetype="dashed", color = "black") + 
  geom_point(aes(x = mat[1,3],y = mat[4,3]), color = "black")

gfg_plot2 <- ggplot(data = df3) +  geom_line(aes(x = Time1,y = RMSE1), color = 2) + 
  geom_point(aes(x = Time1,y = RMSE1, color = "Lanczos"))  + 
  geom_text(aes(x = Time1,y = RMSE1, label = vec_num2),hjust=0.5, vjust=-0.5,color = 2,size = 4) +
  geom_text(aes(x = Time1,y = LogScore1, label = vec_num1),hjust=0.5, vjust=11,color = 2,size = 4) +
  geom_line(aes(x = Time,y = RMSE, color = "Stochastic")) + geom_point(aes(x = Time,y = RMSE), color = 4) +
  geom_text(aes(x = Time,y = RMSE, label = vec_num),hjust=-0.1, vjust=0,color = 4,size = 4) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Time (s)", y = "RMSE") + theme(
    legend.position="bottom",
    legend.title=element_blank(),
    legend.background = element_rect(linetype="solid", 
                                     colour ="black"), text = element_text(size=20))




ggpubr::ggarrange(gfg_plot1,gfg_plot2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")


