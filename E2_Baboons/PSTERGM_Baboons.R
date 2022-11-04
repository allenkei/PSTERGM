library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose()) # choose PSTERGM_Baboons.cpp from the same folder
baboons <- read.delim(file.choose()) # choose RFID_data.txt downloaded from https://osf.io/ufs3y/

# Construct data 
baboons[,4] <- as.Date(baboons[,4], format = "%d/%m/%Y %H:%M")
baboons <- baboons[order(baboons$i,baboons$j,baboons$DateTime),] # 20s interval log
id <- unique(c(baboons$i,baboons$j))
start_date <- as.Date('13/06/2019', format = "%d/%m/%Y")
last_date <- as.Date('10/07/2019', format = "%d/%m/%Y")

y_list <- y_list_all <- list()
list_counter = 1

while(start_date <= last_date){
  
  y <- matrix(0, nrow=length(id), ncol=length(id))
  
  track_i <- track_j <- track_time <- Inf # initialize tracker
  end_date <- start_date + 1
  temp <- baboons[which(baboons[,4] >= start_date & baboons[,4] < end_date),]
  
  for(iter in 1:dim(temp)[1]){
    i <- which(temp[iter,2] == id)
    j <- which(temp[iter,3] == id)
    cur_time <- temp[iter,1]
    
    if(i != track_i | j != track_j){
      track_i <- i
      track_j <- j
      y[i,j] <- y[i,j] + 1
      y[j,i] <- y[j,i] + 1
      track_time <- cur_time
    }else if(i == track_i & j == track_j){
      if(cur_time - track_time > 20){
        y[i,j] <- y[i,j] + 1
        y[j,i] <- y[j,i] + 1
      }
      track_time <- cur_time
    }
    
  }
  
  start_date <- end_date
  # y_list is day 1 to day 23 and y_list_all is day 1 to day 28
  if(list_counter <= 23){y_list[[list_counter]] <- y} 
  y_list_all[[list_counter]] <- y
  
  list_counter <- list_counter + 1
}

rm(temp, y, cur_time, end_date, i, j, iter, last_date, start_date, 
   track_i, track_j, track_time, baboons, list_counter)

node_attr <- rep(0, length(id)) # placeholder
edge_attr <- matrix(0, nrow=length(id),ncol=length(id)) # placeholder

# parameter estimation
n <- length(id)
eta <- rep(0, 8)
eta <- partial_stepping_temporal(20, 100, n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=200)
eta <- newton_raphson_temporal(20, 100, 2*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=200)
eta <- newton_raphson_temporal(10, 1000, 50*n*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=200)
write.csv(eta, 'eta.csv', row.names = F)

# standard error
se_eta <- SE_temporal(1000, 20*n*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=200)
write.csv(se_eta, "se_eta.csv", row.names = F)

# Sampled future networks
y_init <- matrix(1,nrow=n,ncol=n); diag(y_init) <- 0
for(iter in 24:28){
  y_after_holder <- gen_MCMC_valued_CD(100, 200*n*n, y_list_all[[iter-1]], y_init, pi0=0.2, m=200, eta[1:4], eta[5:8], node_attr, edge_attr, iter)
  name <- paste0("y",iter,"_syn.RData")
  save(y_after_holder, file = name) # save the sampled data
}

# load the sampled data
load('y24_syn.RData'); y24_syn <- y_after_holder
load('y25_syn.RData'); y25_syn <- y_after_holder
load('y26_syn.RData'); y26_syn <- y_after_holder
load('y27_syn.RData'); y27_syn <- y_after_holder
load('y28_syn.RData'); y28_syn <- y_after_holder

# calculate the sampled network statistics
y24_syn_features <- gen_features_MCMC_valued(y24_syn, y_list_all, 23, node_attr, edge_attr, length(eta)); rm(y24_syn) 
y25_syn_features <- gen_features_MCMC_valued(y25_syn, y_list_all, 24, node_attr, edge_attr, length(eta)); rm(y25_syn)
y26_syn_features <- gen_features_MCMC_valued(y26_syn, y_list_all, 25, node_attr, edge_attr, length(eta)); rm(y26_syn)
y27_syn_features <- gen_features_MCMC_valued(y27_syn, y_list_all, 26, node_attr, edge_attr, length(eta)); rm(y27_syn) 
y28_syn_features <- gen_features_MCMC_valued(y28_syn, y_list_all, 27, node_attr, edge_attr, length(eta)); rm(y28_syn) 

# Calculate observed network statistics
g_obs <- matrix(0, nrow=5, ncol=length(eta))
for(iter in 24:28){
  aug_dim_holder <- aug_dim_network(y_list_all[[iter-1]], y_list_all[[iter]])
  y_aug <- aug_dim_holder[[1]]; y_dim <- aug_dim_holder[[2]]
  g_obs[iter-23,] <- c(gen_feature_valued(y_aug, node_attr, edge_attr, iter),gen_feature_valued(y_dim, node_attr, edge_attr, iter))
}

rm(aug_dim_holder, y_aug, y_dim, iter)

# generate boxplot for temporal trends
par(mfrow=c(2,4),mar = c(2,2,2,2))
name <- c("Aug: edge sum", "Aug: propensity", "Aug: dispersion", "Aug: transitivity",
          "Dim: edge sum", "Dim: propensity", "Dim: dispersion", "Dim: transitivity")

for(iter in 1:8){
  boxplot(y24_syn_features[,iter],y25_syn_features[,iter],y26_syn_features[,iter],
          y27_syn_features[,iter],y28_syn_features[,iter],main = name[iter], names=24:28)
  lines(1:5,g_obs[,iter], col='red', lwd=2, lty=1)
}

