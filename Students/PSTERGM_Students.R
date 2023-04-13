library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose()) # choose PSTERGM_Students.cpp from the same folder

# data downloaded from http://www.sociopatterns.org/datasets/high-school-contact-and-friendship-networks/
d <- read.csv(file.choose(), sep = ',', header = T)    # High-School_data_2013.csv
g <- read.table(file.choose())                         # metadata_2013.txt
fb <- read.csv(file.choose(), sep = ' ', header = F)   # Facebook-known-pairs_data_2013.csv

# one of the nine classes
class_choice <- "MP"

data <- d[which(d$V4 == class_choice & d$V5 == class_choice),1:3]
facebook <- fb

gender <- g[which(g$V2 == class_choice),c(1,3)]
gender <- gender[-c(29:32),] # remove unknown gender
node_id <- c(gender[,1]) 

colnames(data) <- c("datetime", "sender", "receiver")
colnames(gender) <- c("id", "gender")
colnames(facebook) <- c("id1", "id2", "friend")

##################
# Node attribute #
##################

node_attr <- ifelse(gender[,2] == "M", 0, 1)

################################
# Facebook (undirected binary) #
################################

edge_attr <- matrix(0, nrow=length(node_id), ncol=length(node_id))

for(iter in 1:dim(facebook)[1]){
  track_i <- facebook[iter,1]
  track_j <- facebook[iter,2]
  
  if((track_i %in% node_id) & (track_j %in% node_id)){
    i <- which(track_i == node_id)
    j <- which(track_j == node_id)
    
    if(edge_attr[i,j] == 0 & facebook[iter,3] == 1){
      edge_attr[i,j] <- edge_attr[i,j] + 1
      edge_attr[j,i] <- edge_attr[j,i] + 1
    }
  }
}

#########################################
# Number of Contact (undirected valued) #
#########################################

nrow(data[which(data$sender > data$receiver),]) # should be 0
data <- data[order(data$sender,data$receiver,data$datetime),] # 20s interval log

start_date <- 1385942400 # lubridate::as_datetime(1385942400) # "2013-12-02 UTC"
y_list_all <- list()

for(iter in 1:5){
  
  track_i <- track_j <- track_time <- Inf # initialize tracker
  end_date <- start_date + 24*60*60
  
  temp <- data[which(data$datetime >= start_date & data$datetime < end_date),]
  
  y <- matrix(0, nrow=length(node_id), ncol=length(node_id))
  
  for(row_idx in 1:dim(temp)[1]){
    
    if( temp[row_idx, 2] %in% node_id & temp[row_idx, 3] %in% node_id ){
      i <- which(temp[row_idx, 2] == node_id)
      j <- which(temp[row_idx, 3] == node_id)
      cur_time <- temp[row_idx, 1]
      
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
    
  }
  
  y_list_all[[iter]] <- y
  start_date <- end_date
  
}

rm(temp, y, i,j, iter, row_idx, start_date, end_date, track_i, track_j,
   cur_time, track_time, d, g, fb, class_choice, data, gender, facebook)

######################
# Parameter learning #
######################

n <- length(node_id)
eta_holder <- se_eta_holder <- matrix(0, nrow=4, ncol=12)

for(iter in 2:5){
  y_list <- list()
  y_list[[1]] <- y_list_all[[iter-1]]
  y_list[[2]] <- y_list_all[[iter]]
  eta <- rep(0, 12)
  eta <- partial_stepping_temporal(20, 100, 5*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=0)
  eta <- newton_raphson_temporal(20, 100, 10*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=0)
  eta <- newton_raphson_temporal(10, 1000, 25*n*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=0)
  eta_holder[iter-1,] <- eta
  
  se_eta_holder[iter-1,] <- SE_temporal(1000, 20*n*n, y_list, eta, node_attr, edge_attr, pi0=0.2, m=0)
  
  temp = aug_dim_network(y_list[[1]], y_list[[2]])
  temp = max(temp[[2]]) # calculate m^t for sampling
  
  y_init <- matrix(0,nrow=n,ncol=n) # diag(y_init) <- 0
  y_after_holder <- gen_MCMC_valued_CD(100, 200*n*n, y_list_all[[iter-1]], y_init, pi0=0.2, m=temp, eta[1:6], eta[7:12], node_attr, edge_attr, iter)
  name <- paste0("y",iter,"_syn.RData")
  save(y_after_holder, file = name)
}

write.csv(eta_holder, "eta_holder.csv", row.names = F)
write.csv(se_eta_holder, "se_eta_holder.csv", row.names = F)

############
# Boxplots #
############

load('y2_syn.RData'); y2_syn <- y_after_holder
load('y3_syn.RData'); y3_syn <- y_after_holder
load('y4_syn.RData'); y4_syn <- y_after_holder
load('y5_syn.RData'); y5_syn <- y_after_holder

y2_syn_features <- gen_features_MCMC_valued(y2_syn, y_list_all, 1, node_attr, edge_attr, 12); rm(y2_syn)
y3_syn_features <- gen_features_MCMC_valued(y3_syn, y_list_all, 2, node_attr, edge_attr, 12); rm(y3_syn)
y4_syn_features <- gen_features_MCMC_valued(y4_syn, y_list_all, 3, node_attr, edge_attr, 12); rm(y4_syn) 
y5_syn_features <- gen_features_MCMC_valued(y5_syn, y_list_all, 4, node_attr, edge_attr, 12); rm(y5_syn) 

g_obs <- matrix(0, nrow=4, ncol=12)
for(iter in 2:5){
  aug_dim_holder <- aug_dim_network(y_list_all[[iter-1]], y_list_all[[iter]])
  y_aug <- aug_dim_holder[[1]]
  y_dim <- aug_dim_holder[[2]]
  g_obs[iter-1,] <- c(gen_feature_valued(y_aug, node_attr, edge_attr, iter),gen_feature_valued(y_dim, node_attr, edge_attr, iter))
}

rm(y_aug, y_dim, aug_dim_holder)

name <- c("Inc: edge sum", "Inc: dispersion", "Inc: homophily", "Inc: heterophily", "Inc: facebook", "Inc: transitivity",
          "Dec: edge sum", "Dec: dispersion", "Dec: homophily", "Dec: heterophily", "Dec: facebook", "Dec: transitivity")

# 10 by 4
par(mfrow=c(2,6),mar = c(2,2,2,2))
for(iter in 1:12){
  data <- c(y2_syn_features[,iter],y3_syn_features[,iter],y4_syn_features[,iter],y5_syn_features[,iter],g_obs[,iter])
  boxplot(y2_syn_features[,iter],y3_syn_features[,iter],y4_syn_features[,iter],y5_syn_features[,iter],
          main=name[iter],names = 2:5, ylim=c(min(data),max(data)))
  lines(1:4,g_obs[,iter], col='red', lwd=2)
}

#########################
# Network visualization #
#########################

library(igraph)
library(network)
library(sna)
library(ndtv)

c_scale <- colorRamp(c('gray85','gray65','gray50'))

lb <- 0; ub <- 90
par(mfrow = c(2, 3), mar = c(2,2,2,2))
for(iter in 1:5){
  net2 <- graph_from_adjacency_matrix(y_list_all[[iter]], mode="undirected", weighted = T)
  E(net2)$width <- 8*(E(net2)$weight - lb)/(ub - lb)
  E(net2)$color <- apply( c_scale((E(net2)$weight - lb)/(ub - lb)), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
  V(net2)$color <- ifelse(node_attr==0,'orange','deepskyblue')
  V(net2)$shape <- ifelse(node_attr==0,'square','circle')
  plot(net2, vertex.label=NA, vertex.size=7, vertex.frame.color=NA,
       layout = layout_on_sphere, main = paste('Day',iter));box(lwd=1)
}

net2 <- graph_from_adjacency_matrix(edge_attr, mode="undirected")
E(net2)$width <- 0.3
V(net2)$color <- ifelse(node_attr==0,'orange','deepskyblue')
V(net2)$shape <- ifelse(node_attr==0,'square','circle')
plot(net2, vertex.label=NA, vertex.size=7, vertex.frame.color=NA, 
     layout = layout_on_sphere, main = 'Facebook (binary)');box(lwd=1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("Male", "Female"), 
       col = c("orange","deepskyblue"), border = "black", pch = c(15,16),
       xpd = TRUE, horiz = TRUE, cex = 1, bty = 'n')


