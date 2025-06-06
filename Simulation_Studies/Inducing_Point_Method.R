################################################################################
### Analysis of inducing point methods
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


#######################
### Compare NEGLL
#######################

# Effective ranges
vec_ER <- c(0.5,0.2,0.05)
# Taper Range
taper_range <- 0.016  # 80 average non-zero entries per row
if(Toy){
  taper_range <- 0.052
}
# Inducing point methods
vec_ind_points_method <- c("random","kmeans++","cover_tree")
# Initialize matrices and vectors
mat_NEGLL_FITC <- matrix(0,25,4*3)
colnames(mat_NEGLL_FITC) <- c(paste(vec_ind_points_method,100),
                              paste(vec_ind_points_method,200),
                              paste(vec_ind_points_method,500),
                              paste(vec_ind_points_method,1000))

mat_NEGLL_FSA <- matrix(0,25,4*3)
colnames(mat_NEGLL_FSA) <- c(paste(vec_ind_points_method,100),
                             paste(vec_ind_points_method,200),
                             paste(vec_ind_points_method,500),
                             paste(vec_ind_points_method,1000))

if(Toy){
  mat_NEGLL_FITC <- matrix(0,10,4*3)
  colnames(mat_NEGLL_FITC) <- c(paste(vec_ind_points_method,10),
                                paste(vec_ind_points_method,50),
                                paste(vec_ind_points_method,100),
                                paste(vec_ind_points_method,200))
  
  mat_NEGLL_FSA <- matrix(0,10,4*3)
  colnames(mat_NEGLL_FSA) <- c(paste(vec_ind_points_method,10),
                               paste(vec_ind_points_method,50),
                               paste(vec_ind_points_method,100),
                               paste(vec_ind_points_method,200))
}

l_FITC <- list()
l_FSA <- list()

ntrails <- 25
if(Toy){
  ntrails <- 10
}
ind <- 1
vec_ct_range <- 1/c(12,17,27,40)
if(Toy){
  vec_ct_range <- vec_ct_range*3
}
vec_ind_points <- c(100,200,500,1000)
if(Toy){
  vec_ind_points <- c(10,50,100,200)
}
for (i in vec_ER) {
  
  # Effective Range
  effectiverange <- i
  
  # Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
  arange <- effectiverange/2.7
  
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
  ij <- 0
  for (j in vec_ind_points_method) {
    # Number of inducing points
    for (ji in 1:length(vec_ind_points)) {
      for (ii in 1:ntrails) {
        
        # FITC
        gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,
                            likelihood = likelihood,num_ind_points = vec_ind_points[ji],cover_tree_radius = vec_ct_range[ji],                 
                            gp_approx = "fitc",ind_points_selection = j,
                            matrix_inversion_method = "cholesky",seed = ii*10)
        
        mat_NEGLL_FITC[ii,ij+ji] <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
        if(!Toy){
        # FSA
          gp_model <- GPModel(gp_coords = coords_train, cov_function = "matern",cov_fct_shape = 1.5,cov_fct_taper_shape = 2,
                              likelihood = likelihood,num_ind_points = vec_ind_points[ji],cov_fct_taper_range = taper_range,                          
                              gp_approx = "full_scale_tapering", cover_tree_radius = vec_ct_range[ji],ind_points_selection = j,
                              matrix_inversion_method = "cholesky",seed = ii*10)
        
          mat_NEGLL_FSA[ii,ij+ji] <- gp_model$neg_log_likelihood(y = y_train,cov_pars = init_cov_pars)
        }
      } 
    }
    ij <- ij + length(vec_ind_points)
    
  }
  l_FITC[[ind]] <- mat_NEGLL_FITC
  l_FSA[[ind]] <- mat_NEGLL_FSA
  ind <- ind+1
  
}


###################
### Plots
###################

# FITC
FITC_plots <- vector('list', 3)
for (i in 1:3) {
  
  mat_NEGLL <- l_FITC[[i]]
  if(Toy){
    mat_NEGLL2 <- cbind(rep(c(rep(10,10),rep(50,10),rep(100,10),rep(200,10)),3),c(rep("Random",4*10),rep("kMeans++",4*10),rep("CoverTree",4*10)),
                        as.numeric(c(as.numeric(mat_NEGLL[1:10,1]),as.numeric(mat_NEGLL[1:10,2]),as.numeric(mat_NEGLL[1:10,3]),as.numeric(mat_NEGLL[1:10,4]),
                                     as.numeric(mat_NEGLL[1:10,5]),as.numeric(mat_NEGLL[1:10,6]),as.numeric(mat_NEGLL[1:10,7]),as.numeric(mat_NEGLL[1:10,8]),
                                     as.numeric(mat_NEGLL[1:10,9]),as.numeric(mat_NEGLL[1:10,10]),as.numeric(mat_NEGLL[1:10,11]),as.numeric(mat_NEGLL[1:10,12]))))
  } else {
    mat_NEGLL2 <- cbind(rep(c(rep(100,25),rep(200,25),rep(500,25),rep(1000,25)),3),c(rep("Random",4*25),rep("kMeans++",4*25),rep("CoverTree",4*25)),
                        as.numeric(c(as.numeric(mat_NEGLL[1:25,1]),as.numeric(mat_NEGLL[1:25,2]),as.numeric(mat_NEGLL[1:25,3]),as.numeric(mat_NEGLL[1:25,4]),
                                     as.numeric(mat_NEGLL[1:25,5]),as.numeric(mat_NEGLL[1:25,6]),as.numeric(mat_NEGLL[1:25,7]),as.numeric(mat_NEGLL[1:25,8]),
                                     as.numeric(mat_NEGLL[1:25,9]),as.numeric(mat_NEGLL[1:25,10]),as.numeric(mat_NEGLL[1:25,11]),as.numeric(mat_NEGLL[1:25,12]))))
  }
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat) <- c("x","Method","value")
  data_mat$value <- as.numeric(data_mat$value)
  
  FITC_plots[[i]] <- local({
    ggplot(data_mat,            
           aes(x = reorder(x,-value),
               y = value,
               fill = Method)) +  geom_boxplot() +  scale_colour_manual(values = c(rgb(1,0,0),
                                                                                   rgb(0,0,1),
                                                                                   rgb(0,1,0)),name = "Method", 
                                                                        labels = c("Random","kMeans++","CoverTree"))+
      theme(
        legend.position="top",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = "", y = "NEGLL")
  })
  
  
}
ggarrange( FITC_plots[[1]],  FITC_plots[[2]],  FITC_plots[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

if(!Toy){
# FSA
FSA_plots <- vector('list', 3)
for (i in 1:3) {
  
  mat_NEGLL <- l_FSA[[i]]
  if(Toy){
    mat_NEGLL2 <- cbind(rep(c(rep(10,10),rep(50,10),rep(100,10),rep(200,10)),3),c(rep("Random",4*10),rep("kMeans++",4*10),rep("CoverTree",4*10)),
                        as.numeric(c(as.numeric(mat_NEGLL[1:10,1]),as.numeric(mat_NEGLL[1:10,2]),as.numeric(mat_NEGLL[1:10,3]),as.numeric(mat_NEGLL[1:10,4]),
                                     as.numeric(mat_NEGLL[1:10,5]),as.numeric(mat_NEGLL[1:10,6]),as.numeric(mat_NEGLL[1:10,7]),as.numeric(mat_NEGLL[1:10,8]),
                                     as.numeric(mat_NEGLL[1:10,9]),as.numeric(mat_NEGLL[1:10,10]),as.numeric(mat_NEGLL[1:10,11]),as.numeric(mat_NEGLL[1:10,12]))))
  } else {
    mat_NEGLL2 <- cbind(rep(c(rep(100,25),rep(200,25),rep(500,25),rep(1000,25)),3),c(rep("Random",4*25),rep("kMeans++",4*25),rep("CoverTree",4*25)),
                        as.numeric(c(as.numeric(mat_NEGLL[1:25,1]),as.numeric(mat_NEGLL[1:25,2]),as.numeric(mat_NEGLL[1:25,3]),as.numeric(mat_NEGLL[1:25,4]),
                                     as.numeric(mat_NEGLL[1:25,5]),as.numeric(mat_NEGLL[1:25,6]),as.numeric(mat_NEGLL[1:25,7]),as.numeric(mat_NEGLL[1:25,8]),
                                     as.numeric(mat_NEGLL[1:25,9]),as.numeric(mat_NEGLL[1:25,10]),as.numeric(mat_NEGLL[1:25,11]),as.numeric(mat_NEGLL[1:25,12]))))
  }
  data_mat <- as.data.frame(mat_NEGLL2)
  rownames(data_mat) <- NULL
  colnames(data_mat) <- c("x","Method","value")
  data_mat$value <- as.numeric(data_mat$value)
  
  FSA_plots[[i]] <- local({
    ggplot(data_mat,            
           aes(x = reorder(x,-value),
               y = value,
               fill = Method)) +  geom_boxplot() +  scale_colour_manual(values = c(rgb(1,0,0),
                                                                                   rgb(0,0,1),
                                                                                   rgb(0,1,0)),name = "Method", 
                                                                        labels = c("Random","kMeans++","CoverTree"))+
      theme(
        legend.position="top",
        legend.background = element_rect(linetype="solid", 
                                         colour ="black")
      ) +
      labs(x = "", y = "NEGLL")
  })
  
  
}
ggarrange( FSA_plots[[1]],  FSA_plots[[2]],  FSA_plots[[3]], ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
}
