################################################################################
### Real world Data Analysis
################################################################################

#######
## Packages
#######

# Package names
packages <- c("dplyr","ggpubr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

library(gpboost)

# Function for Simulating Data
source("https://raw.github.com/TimGyger/iterativeFSA/master/Data/Load_Data.R")

#####################################################
# Comparison
#####################################################

dat_train <- as.matrix(data_train)
dat_test <- as.matrix(data_test)
y_test <- dat_test[,3]
mat <- matrix(0,11,4)
colnames(mat) <- c("Iterative FSA","Exact FSA","FITC","Tapering")
rownames(mat) <- c("beta_0","beta_1","beta_2","sigma","sigma_1","rho","time1","RMSE","Log-Score","CRPS","time2")
###########################
### Iterative FSA
###########################

a <- Sys.time()
gp_model <- fitGPModel(X = as.matrix(cbind(rep(1,nrow(dat_train)),dat_train[,1:2])),gp_coords = dat_train[,1:2], cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "iterative",
                       likelihood = "gaussian",seed = 10, num_ind_points = 500,cov_fct_taper_range = 7500,cov_fct_taper_shape = 2, gp_approx = "full_scale_tapering",
                       y = dat_train[,3],params = list(maxit=1000,trace=TRUE,cg_delta_conv = 1,cg_preconditioner_type = "predictive_process_plus_diagonal",
                                                                                                  optimizer_cov = "lbfgs",
                                                                                                  cg_max_num_it = 1000,cg_max_num_it_tridiag = 1000,num_rand_vec_trace = 50,
                                                                                                  seed_rand_vec_trace = 10,reuse_rand_vec_trace = T))


mat[1:7,1] <- c(gp_model$get_coef(),gp_model$get_cov_pars(),Sys.time()-a)

aa <- Sys.time()
gp_model$set_prediction_data(cg_delta_conv_pred = 1, nsim_var_pred = 500)
pred <- predict(X_pred = as.matrix(cbind(rep(1,nrow(dat_test)),dat_test[,1:2])),gp_model, gp_coords_pred = as.matrix(dat_test[,1:2]), predict_var = T,y = dat_train[,3])
time2 <- Sys.time()-aa

mat[8,1] <- sqrt(mean((pred$mu-y_test)^2))
mat[9,1] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
mat[10,1] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
mat[11,1] <- time2

pred1 <- pred

###########################
### Exact FSA
###########################

a <- Sys.time()
gp_model <- fitGPModel(X = as.matrix(cbind(rep(1,nrow(dat_train)),dat_train[,1:2])),gp_coords = dat_train[,1:2], cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = "gaussian",seed = 10, num_ind_points = 500,cov_fct_taper_range = 7500,cov_fct_taper_shape = 2, gp_approx = "full_scale_tapering",
                       y = dat_train[,3],params = list(maxit=1000,trace=TRUE,optimizer_cov = "lbfgs"))


mat[1:7,2] <- c(gp_model$get_coef(),gp_model$get_cov_pars(),Sys.time()-a)

aa <- Sys.time()
pred <- predict(X_pred = as.matrix(cbind(rep(1,nrow(dat_test)),dat_test[,1:2])),gp_model, gp_coords_pred = as.matrix(dat_test[,1:2]), predict_var = T,y = dat_train[,3])
time2 <- Sys.time()-aa

mat[8,2] <- sqrt(mean((pred$mu-y_test)^2))
mat[9,2] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
mat[10,2] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
mat[11,2] <- time2

pred2 <- pred

###########################
### FITC
###########################

a <- Sys.time()
gp_model <- fitGPModel(X = as.matrix(cbind(rep(1,nrow(dat_train)),dat_train[,1:2])),gp_coords = dat_train[,1:2], cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = "gaussian",seed = 10, num_ind_points = 500,gp_approx = "FITC",
                       y = dat_train[,3],params = list(maxit=1000,trace=TRUE, optimizer_cov = "lbfgs"))

mat[1:7,3] <- c(gp_model$get_coef(),gp_model$get_cov_pars(),Sys.time()-a)

aa <- Sys.time()
pred <- predict(X_pred = as.matrix(cbind(rep(1,nrow(dat_test)),dat_test[,1:2])),gp_model, gp_coords_pred = as.matrix(dat_test[,1:2]), predict_var = T,y = dat_train[,3])
time2 <- Sys.time()-aa

mat[8,3] <- sqrt(mean((pred$mu-y_test)^2))
mat[9,3] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
mat[10,3] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
mat[11,3] <- time2

pred3 <- pred


###########################
### Tapering
###########################

a <- Sys.time()
gp_model <- fitGPModel(X = as.matrix(cbind(rep(1,nrow(dat_train)),dat_train[,1:2])),gp_coords = dat_train[,1:2], cov_function = "matern",cov_fct_shape = 1.5,matrix_inversion_method = "cholesky",
                       likelihood = "gaussian",seed = 10, cov_fct_taper_range = 7500,cov_fct_taper_shape = 2, gp_approx = "tapering",
                       y = dat_train[,3],params = list(maxit=1000,trace=TRUE,optimizer_cov = "gradient_descent",optimizer_coef = "lbfgs",init_cov_pars = c(21.3375,21.3375,207531)))

mat[1:7,4] <- c(gp_model$get_coef(),gp_model$get_cov_pars(),Sys.time()-a)

aa <- Sys.time()
pred <- predict(X_pred = as.matrix(cbind(rep(1,nrow(dat_test)),dat_test[,1:2])),gp_model, gp_coords_pred = as.matrix(dat_test[,1:2]), predict_var = T,y = dat_train[,3])
time2 <- Sys.time()-aa

mat[8,4] <- sqrt(mean((pred$mu-y_test)^2))
mat[9,4] <- -mean(dnorm(y_test,pred$mu,sqrt(pred$var), log = T))
help_vec <- dnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
help_vec2 <- pnorm((y_test-pred$mu)/sqrt(pred$var),rep(0,length(pred$mu)),rep(1,length(pred$mu)))
mat[10,4] <- -mean(sqrt(pred$var)*(1/sqrt(pi)-2*help_vec-((y_test-pred$mu)/sqrt(pred$var))*(2*help_vec2-1)))
mat[11,4] <- time2

pred4 <- pred

########################
### Table
########################

print(mat)

#########################
### Plots
#########################

par(mfrow = c(2,2),mai = c(0.2, 0.25, 0.25, 0.01))

quilt.plot(dat_test[,1],dat_test[,2],y_test,nx = 200,ny = 200,legend.width = 0.5,xlab = "",ylab = "Northing (m)",main="Iterative-FSA",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred1$mu,nx = 200,ny = 200,legend.width = 0.5,xlab = "",ylab = "",main="Cholesky-FSA",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred3$mu,nx = 200,ny = 200,legend.width = 0.5,xlab = "Easting (m)",ylab = "Northing (m)",main="Tapering",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred4$mu,nx = 200,ny = 200,legend.width = 0.5,xlab = "Easting (m)",ylab = "",main="FITC",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)

quilt.plot(dat_test[,1],dat_test[,2],pred1$var,nx = 200,ny = 200,legend.width = 0.5,xlab = "",ylab = "Northing (m)",main="Iterative-FSA",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred2$var,nx = 200,ny = 200,legend.width = 0.5,xlab = "",ylab = "",main="Cholesky-FSA",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred3$var,nx = 200,ny = 200,legend.width = 0.5,xlab = "Easting (m)",ylab = "Northing (m)",main="Tapering",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)
quilt.plot(dat_test[,1],dat_test[,2],pred4$var,nx = 200,ny = 200,legend.width = 0.5,xlab = "Easting (m)",ylab = "",main="FITC",cex.axis=1.8,axis.args=list(cex.axis=2),cex.lab=2, cex.main=2)