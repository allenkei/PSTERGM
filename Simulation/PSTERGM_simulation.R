library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose()) # choose PSTERGM_S2.cpp from the same folder

# specify T = 5, 10, 15 for this experiment
length_of_T <- 10
n <- 50

# network50.csv has 50 nodes and network100.csv has 100 nodes
network0 <- as.matrix(read.csv(file.choose())) # choose network50.csv or network100.csv from the same folder

node_attr <- rep(1,n)
eta_atleast <- c(-2,2,1,1) # eta_plus
eta_atmost <- c(-1,2,1,1) # eta_minus


par_iter <- 50
eta_holder <- se_holder <- matrix(0, nrow=par_iter, ncol=length(c(eta_atleast,eta_atmost)))

for(i in 41:50){
  y_list <- list(); y_list[[1]] <- network0
  y_init <- matrix(0,nrow=n,ncol=n);
  
  # generate sampled data
  for(t in 2:length_of_T){
    MCMC_list <- gen_MCMC_valued_CD(1, 20*n*n, y_list[[t-1]], y_init, pi0=0.2, m=3, eta_atleast, eta_atmost, node_attr)
    y_list[[t]] <- MCMC_list[[1]]
  }
  
  save(y_list, file = paste0("T",length_of_T,"_n",n,"_iter",i,".Rdata")) 
  
  # parameter estimation
  eta <- rep(0, length(c(eta_atleast,eta_atmost)))
  eta <- partial_stepping_temporal(20, 100, n, y_list, eta, node_attr, pi0=0.2, m=3)
  eta <- newton_raphson_temporal(5, 1000, n, y_list, eta, node_attr, pi0=0.2, m=3)
  eta_holder[i,] <- eta
  se_holder[i,] <- SE_temporal(100, 1000, y_list, eta, node_attr, pi0=0.2, m=3) # use learned eta
}


write.csv(eta_holder, 'eta_holder5.csv', row.names = F)
write.csv(se_holder, 'se_holder5.csv', row.names = F)


rm(y_init, eta, MCMC_list, y_list, i, node_attr, t, network0)


# Result for Boxplot

y_temp2 <- y_temp3 <- y_temp4 <- y_temp5 <- y_temp6 <- y_temp7 <- y_temp8 <- y_temp9 <- y_temp10 <- matrix(0,nrow=par_iter,ncol=4)

for(i in 1:par_iter){
  load(paste0("T",length_of_T,"_n",n,"_iter",i,".Rdata"))
  
  y_temp2[i,] <- c(gen_feature_valued(y_list[[2]], node_attr))
  y_temp3[i,] <- c(gen_feature_valued(y_list[[3]], node_attr))
  y_temp4[i,] <- c(gen_feature_valued(y_list[[4]], node_attr))
  y_temp5[i,] <- c(gen_feature_valued(y_list[[5]], node_attr))
  y_temp6[i,] <- c(gen_feature_valued(y_list[[6]], node_attr))
  y_temp7[i,] <- c(gen_feature_valued(y_list[[7]], node_attr))
  y_temp8[i,] <- c(gen_feature_valued(y_list[[8]], node_attr))
  y_temp9[i,] <- c(gen_feature_valued(y_list[[9]], node_attr))
  y_temp10[i,] <- c(gen_feature_valued(y_list[[10]], node_attr))
}


par(mfrow=c(1,4),mar = c(2,2,2,2))
name <- c("Edge sum", "Zeros", "Mutuality", "Transitivity")

for(iter in 1:4){
  boxplot(y_temp2[,iter],y_temp3[,iter],y_temp4[,iter],y_temp5[,iter],
          y_temp6[,iter],y_temp7[,iter],y_temp8[,iter],y_temp9[,iter],y_temp10[,iter],
          main = name[iter], names=2:10)
}



# Result for Table

round(abs(apply(eta_holder,2,median) - c(eta_atleast,eta_atmost)),4)
round(apply(eta_holder,2,sd),3)


# Result for Plot

eta_true <- c(eta_atleast, eta_atmost)
CI_data <- data.frame(matrix(nrow = par_iter*8, ncol = 6))
colnames(CI_data) <- c('ind', 'est', 'lower', 'upper', 'true', 'col')
name <- c("Inc: edge sum", "Inc: zeros", "Inc: mutuality", "Inc: transitivity",
          "Dec: edge sum", "Dec: zeros", "Dec: mutuality", "Dec: transitivity")



CI_data$ind <- rep(1:par_iter,8)
CI_data$true <- rep(eta_true, each=par_iter)


for(i in 1:8){
  CI_data[ ((i-1)*par_iter+1):(i*par_iter), 2] <- eta_holder[,i]
  CI_data[ ((i-1)*par_iter+1):(i*par_iter), 3] <- eta_holder[,i] - 1.96 * se_holder[,i]
  CI_data[ ((i-1)*par_iter+1):(i*par_iter), 4] <- eta_holder[,i] + 1.96 * se_holder[,i]
}


CI_data$col <- ifelse( (CI_data$lower<CI_data$true)&(CI_data$upper>CI_data$true), 'black', 'red')



library(ggplot2)
library(cowplot)
p <- list()
for(i in 1:8){
  temp <- CI_data[((i-1)*par_iter+1):(i*par_iter), ]
  p[[i]] <- ggplot(temp, aes(ind, est, color=col)) +        # ggplot2 plot with confidence intervals
    geom_point(size=0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),linewidth=0.2) +
    geom_hline(yintercept=eta_true[i], color = "blue") +
    scale_colour_identity() + theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5)) +
    ggtitle(name[i])
}



plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]],
          ncol = 4, align = "v")




