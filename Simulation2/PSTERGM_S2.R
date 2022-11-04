library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose()) # choose PSTERGM_S2.cpp from the same folder

# specify T = 5, 10, 15 for this experiment
length_of_T <- 5 

# specify n = 50, 100 for this experiment
# network50.csv has 50 nodes and network100.csv has 100 nodes
network_before <- as.matrix(read.csv(file.choose())) # choose network50.csv or network100.csv from the same folder

n <- dim(network_before)[1]; node_attr <- rep(1,n)
eta_atleast <- c(-2,1,1) # eta_plus
eta_atmost <- c(-1,1,1) # eta_minus

y_list <- list(); y_list[[1]] <- network_before
y_init <- matrix(1,nrow=n,ncol=n); diag(y_init) <- 0;

# generate sampled data
for(t in 2:length_of_T){
  MCMC_list <- gen_MCMC_valued_CD(1, 50*n*n, network_before, y_init, pi0=0.2, m=3, eta_atleast, eta_atmost, node_attr)
  y_list[[t]] <- network_before <- MCMC_list[[1]]
}

#save(y_list, file = "T5_node50.Rdata") # change file name based on specification
#load("/home/yikkie/Paper/E2/T5_node50/T5_node50.Rdata")

# parameter estimation
par_iter <- 100
eta_holder <- matrix(0,nrow=par_iter,ncol=6)
for(i in 1:par_iter){
  
  eta <- c(0, 0, 0, 0, 0, 0)
  eta <- partial_stepping_temporal(20, 100, n, y_list, eta, node_attr, pi0=0.2, m=3)
  eta_holder[i,] <- newton_raphson_temporal(5, 1000, n, y_list, eta, node_attr, pi0=0.2, m=3)
  
}

round(abs(apply(eta_holder,2,median) - c(eta_atleast,eta_atmost)),3)
round(apply(eta_holder,2,sd),3)

#write.csv(eta_holder, 'eta_holder.csv', row.names = F)

