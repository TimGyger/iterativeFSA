################################################################################
### Preconditioner Comparison (Vecchia)
################################################################################


#######
## Packages
#######

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../Packages.R")

################################################################################
#Generate data

make_data <- function(n, 
                      likelihood = "bernoulli_logit",
                      sigma2=1, #marginal variance of GP
                      rho=0.1, #range parameter
                      cov_function="exponential",
                      matern_nu=2.5,
                      with_fixed_effects=FALSE,
                      beta,
                      gamma_shape = 1,
                      p) {
  
  #Simulate spatial Gaussian process
  coords <- matrix(runif(n*2),ncol=2)
  if (cov_function == "exponential"){
    RFmodel <- RMexp(var=sigma2, scale=rho)
  } else if(cov_function == "matern") {
    RFmodel <- RMmatern(matern_nu, notinvnu=TRUE, var=sigma2, scale=rho)
  } else {
    stop("cov_function not supported")
  }
  sim <- RFsimulate(RFmodel, x=coords)
  eps <- sim$variable1
  eps <- eps - mean(eps)
  if(with_fixed_effects){
    #   Simulate fixed effects
    X <- cbind(rep(1,n),matrix(runif(n*p)-0.5,ncol=p))
    f <- X %*% beta    
  } else{
    X <- NA
    f <- 0
  }
  b <- f+eps
  #Simulate response variable
  if (likelihood == "bernoulli_probit") {
    probs <- pnorm(b)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(b)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(b)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    y <- qgamma(runif(n), rate = gamma_shape * exp(-b), shape = gamma_shape) #E[y] = exp(b)
  } else if (likelihood == "gaussian") {
    mu <- b
    y <- rnorm(n,sd=0.05) + mu
  }
  list(y=y, coords=coords, b=b, X=X)
}

### Toy example
Toy <- TRUE

sigma2_true <- 1
rho_true <- 0.25 #other ranges considered: 0.01, 0.25
true_covpars <- c(sigma2_true, rho_true)
n <- 100000
if (Toy){
  n <- 10000
}
NN <- 20

set.seed(1)
mydata <- make_data(n=n,
                    likelihood = "bernoulli_logit",
                    sigma2=sigma2_true,
                    rho=rho_true,
                    cov_function="matern",
                    matern_nu=1.5,
                    with_fixed_effects=FALSE)

################################################################################
NUM_RAND_VEC_TRACE <- c(10, 20, 50, 100)
PRECONDITIONER <- c("piv_chol_on_Sigma", "Sigma_inv_plus_BtWB","predictive_process_plus_diagonal_200",
                    "predictive_process_plus_diagonal_500", "vecchia_response", "incomplete_cholesky")
n_rep <- 10

VLresult <- NA
VLtime <- NA

Itresults <- data.frame(matrix(nrow=length(NUM_RAND_VEC_TRACE)*length(PRECONDITIONER)*n_rep,ncol = 4))
colnames(Itresults) <- c("preconditioner", "t", "negLL", "time")

################################################################################
#Iterative-VL

i = 1
for(p in 1:length(PRECONDITIONER)){
  for(t in 1:length(NUM_RAND_VEC_TRACE)){
    for(r in 1:n_rep){
      print(p)
      print(t)
      print(r)
      Itmodel <- GPModel(gp_coords = mydata$coords[,1:2],
                         cov_function = "matern",
                         cov_fct_shape=1.5,
                         likelihood="bernoulli_logit",
                         matrix_inversion_method = "iterative",
                         gp_approx = "vecchia",
                         vecchia_ordering = "random",
                         num_neighbors=NN)
      if (p == 3){
        piv_chol_rank <- 200
        if(Toy){
          piv_chol_rank <- 50
        }
        Itmodel$set_optim_params(params = list(maxit=1,
                                               trace=TRUE,
                                               num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                               cg_preconditioner_type = "predictive_process_plus_diagonal",
                                               seed_rand_vec_trace = i,fitc_piv_chol_preconditioner_rank = piv_chol_rank))
      } else if (p == 4){
        Itmodel$set_optim_params(params = list(maxit=1,
                                               trace=TRUE,
                                               num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                               cg_preconditioner_type = "predictive_process_plus_diagonal",
                                               seed_rand_vec_trace = i,fitc_piv_chol_preconditioner_rank = 500))
      } else {
        Itmodel$set_optim_params(params = list(maxit=1,
                                               trace=TRUE,
                                               num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                               cg_preconditioner_type = PRECONDITIONER[p],
                                               seed_rand_vec_trace = i))
      }
      
      
      Itresults$preconditioner[i] <- PRECONDITIONER[p]
      Itresults$t[i] <- NUM_RAND_VEC_TRACE[t]
      Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=true_covpars, y=mydata$y))[3]
      
      i = i+1
      gc()
    }
  }
}


################################################################################
#Cholesky-VL

VLmodel <- GPModel(gp_coords = mydata$coords[,1:2],
                   cov_function = "matern",
                   cov_fct_shape=1.5,
                   likelihood="bernoulli_logit",
                   matrix_inversion_method = "cholesky",
                   gp_approx = "vecchia",
                   vecchia_ordering = "random",
                   num_neighbors=NN)

VLmodel$set_optim_params(params = list(maxit=1,
                                       trace=TRUE))

VLtime <- system.time(VLresult <- VLmodel$neg_log_likelihood(cov_pars=true_covpars, y=mydata$y))[3]

################################################################################
# Plotting
################################################################################

library(ggplot2)
library(grid)

#Renaming
Itresults$preconditioner[Itresults$preconditioner=="Sigma_inv_plus_BtWB"] <- "P[VADU]"
Itresults$preconditioner[Itresults$preconditioner=="piv_chol_on_Sigma"] <- "P[LRAC]"
Itresults$preconditioner[Itresults$preconditioner=="predictive_process_plus_diagonal_200"] <- "P[FITC_200]"
Itresults$preconditioner[Itresults$preconditioner=="predictive_process_plus_diagonal_500"] <- "P[FITC_500]"
Itresults$preconditioner[Itresults$preconditioner=="vecchia_response"] <- "P[Vecchia]"
Itresults$preconditioner[Itresults$preconditioner=="incomplete_cholesky"] <- "P[Inc_Chol]"
Itresults$preconditioner <- factor(Itresults$preconditioner, levels = c("P[VADU]", "P[LRAC]","P[FITC_200]","P[FITC_500]",
                                                                        "P[Vecchia]","P[Inc_Chol]"), ordered=TRUE)
Itresults$t <- as.factor(Itresults$t)

p1 <- ggplot(Itresults, aes(x=t, y=negLL, fill=preconditioner)) + 
  geom_hline(yintercept=VLresult, linetype = "dashed") +  
  geom_boxplot() + labs(fill  = "") + ylab("log-likelihood") +
  scale_fill_brewer(type = "qual", palette=6, labels = c("VADU", "Pivoted Cholesky", "FITC","FITC_500","Observable Vecchia","Incomplete Cholesky")) +
  theme_bw() + theme(axis.title.x=element_blank(),
                     text = element_text(size=20),
                     axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank(), 
                     legend.position = "top", 
                     axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + scale_y_continuous()

p2 <- ggplot(Itresults, aes(x=t, y=time, color=preconditioner, shape=preconditioner)) +
  stat_summary(aes(group = preconditioner), fun = mean, geom = 'line', size=1) + 
  stat_summary(aes(group = preconditioner), fun = mean, geom = 'point', size=2) +
  scale_color_brewer(type = "qual", palette=6) +labs(color = "") + 
  scale_shape_manual(values = c(1,2,3,4,5,6), labels = scales::parse_format()) +
  theme_bw() + theme(legend.position = "none",text = element_text(size=20), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  ylab("Time (s)") +
  xlab("Number of sample vectors") 

if (!Toy){
  p2 <- p2 + scale_y_log10(breaks = c(2,5,10,15,20,30,50,100)) +
    annotate(geom="text", x=0.8, y=95, label=paste("Cholesky:",round(VLtime),"s"),
             color="black",size = 5) 
}

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), 
                ggplotGrob(p2), size = "first"))
