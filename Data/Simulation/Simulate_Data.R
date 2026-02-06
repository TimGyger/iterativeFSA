#######################################################
### Simulate Data
#######################################################

### Function for simulating Data
sim_data <- function(n,smoothness = 3/2, Covfunct = "matern",range_param,sigma_param,sigma_error,seed) {
  #   Simulate spatial Gaussian process
  sigma2_1 <- sigma_param # marginal variance of GP
  rho <- range_param # range parameter
  coords <- matrix(runif(n*2),ncol=2)
  set.seed(seed)
  if (Covfunct == "matern"){
    RFmodel <- RMmatern(smoothness, var = sigma2_1, scale = rho)
  } else {
    RFmodel <- RMexp(var=sigma2_1, scale=rho) 
  }
  sim <- RFsimulate(RFmodel, x=coords)
  eps <- sim$variable1
  eps <- eps - mean(eps)
  #   Simulate response variable
  mu <- eps
  y <- rnorm(n,sd=sqrt(sigma_error)) + mu
  list(y=y, coords = coords)
}

### Example
example <- F
if (example){
  vec_ranges <- c(0.2,0.05,0.01)
  par(mfrow = c(1,3))
  for (i in 1:3) {
    ## Covariance Function
    # sigma
    sigma <-  1
    # Covariance Function
    Covfunct <-  "Matern"
    # Smoothness Parameter
    smoothness <-  3/2
    # Effective Range
    effectiverange <- vec_ranges[i]
    # Range Parameter (see Table 1 in Jointly Specified Spatial Priors for Bayesian Models of Crash Frequency)
    if (Covfunct == "Exponential"){
      arange <- effectiverange/3
    } else if(Covfunct == "Matern"){
      arange <- effectiverange/4.7439
    }
    print(arange)
    # Variance of error term
    sigma_error = 1
    
    ## Data
    # number of data points
    n <- 100000
    # Ordering Data points (distance to (0,0))
    order <- T
    set.seed(1)
    simdata <- sim_data(n = n,smoothness = smoothness, Covfunct = Covfunct,range_param = arange,sigma_param = sigma, sigma_error = sigma_error, seed = 1)
    
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
    legend_show <- F
    if(i == 3){
      legend_show <- T
    }
    quilt.plot(X[,1],X[,2],Y,nx = 200,ny = 200,add.legend = legend_show,cex.axis=1.8,axis.args=list(cex.axis=2))
    grid()
  }
}
