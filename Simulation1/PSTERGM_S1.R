library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose()) # choose PSTERGM_S1.cpp from the same folder
network_before <- as.matrix(read.csv(file.choose())) # choose network_before.csv from the same folder

n <- 50; node_attr <- rep(1,n)
y_init <- matrix(1,nrow=n,ncol=n); diag(y_init) <- 0
eta_atleast <- c(-2,1,1) # eta_plus
eta_atmost <- c(-1,1,1) # eta_minus

# generate sampled data
MCMC_list <- gen_MCMC_valued_CD(100, 10*n*n, network_before, y_init, pi0=0.2, m=3, eta_atleast, eta_atmost, node_attr)

#save(MCMC_list, file = "100_y2_50nodes.RData")
#load("100_y2_50nodes.RData")

# parameter estimation
eta_holder <- matrix(0,nrow=100,ncol=6)
for(iter in 1:100){

  y_list <- list(network_before, MCMC_list[[iter]])
  
  eta <- c(0, 0, 0, 0, 0, 0)
  eta <- partial_stepping_temporal(20, 100, n, y_list, eta, node_attr, pi0=0.2, m=3)
  eta_holder[iter,] <- newton_raphson_temporal(5, 1000, n, y_list, eta, node_attr, pi0=0.2, m=3)
  
}

round(abs(apply(eta_holder,2,median) - c(eta_atleast, eta_atmost)),3)
round(apply(eta_holder,2,sd),3)

#write.csv(eta_holder, "eta_holder.csv", row.names = F)
