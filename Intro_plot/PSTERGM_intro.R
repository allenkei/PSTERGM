obs <- c(2,5,7,4,3,5,9,7,5,5)

plus <- minu <- c()
for(i in 2:length(obs)){
  plus <- c(plus, max(obs[i],obs[i-1]))
  minu <- c(minu, min(obs[i],obs[i-1]))
}

x <- seq(1.5, length(obs)-0.5, by=1)

par(mar=c(4.1, 4.1, 0.5, 0.5), xpd=TRUE)
plot(1:length(obs), obs, type='l', xlab='time point', ylab = 'value (dyad)', 
     lwd = 2, xaxt='n', yaxt='n',pch=19, cex = 0.5)


xtick<-1:length(obs)
axis(side=1, at=xtick)
axis(side=2, at=c(2,4,6,8))
lines(x, plus, col = 'red', type='b', pch=19, lty=2)
lines(x, minu, col = 'blue', type='b', pch=19, lty=2)

lines(c(2.5,3.5), c(7,7), col = 'red', lwd = 3)
lines(c(6.5,7.5), c(9,9), col = 'red', lwd = 3)
lines(c(4.5,5.5), c(3,3), col = 'blue', lwd = 3)


for(i in c(1,3,5,7,9)){
  rect(xleft = i, xright = i+1, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("grey", alpha = 0.3))
}




library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp(file.choose())
baboons <- read.delim("https://raw.githubusercontent.com/allenkei/Network_Data/main/Baboons/RFID_data.txt") # choose RFID_data.txt downloaded from https://osf.io/ufs3y/

# Construct data 
baboons[,4] <- as.Date(baboons[,4], format = "%d/%m/%Y %H:%M")
baboons <- baboons[order(baboons$i,baboons$j,baboons$DateTime),] # 20s interval log
id <- unique(c(baboons$i,baboons$j))
start_date <- as.Date('13/06/2019', format = "%d/%m/%Y")
last_date <- as.Date('10/07/2019', format = "%d/%m/%Y")

y_list_all <- list()
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
  y_list_all[[list_counter]] <- y
  list_counter <- list_counter + 1
}


node_attr <- rep(0, length(id)) # placeholder
edge_attr <- matrix(0, nrow=length(id),ncol=length(id)) # placeholder
rm(temp, y, cur_time, end_date, i, j, iter, last_date, start_date, 
   track_i, track_j, track_time, baboons, list_counter, id)



holder <- y_aug_holder <- y_dim_holder <- c()
for(i in 23:28){ holder <- c(holder, sum(y_list_all[[i]])/2) }

for(i in 24:28){
  temp <- aug_dim_network(y_list_all[[i-1]], y_list_all[[i]])
  y_aug_holder <- c(y_aug_holder, sum(temp[[1]])/2)
  y_dim_holder <- c(y_dim_holder, sum(temp[[2]])/2)
}


par(mar=c(4, 4, 0.5, 0.5), xpd=TRUE)

plot(23:28, holder, type='l', ylab = 'edge sum (network)', xlab = "day", lwd = 2, yaxt='n',
     ylim=c(min(holder, y_aug_holder, y_dim_holder),max(holder, y_aug_holder, y_dim_holder)))
axis(side=2, at=c(700,900,1100,1300))
lines(23.5:27.5, y_aug_holder, col = 'red', type='b', lty=2, pch=19)
lines(23.5:27.5, y_dim_holder, col = 'blue', type='b', lty=2, pch=19)
points(25, holder[3], pch=19)
points(26, holder[4], pch=19)
segments(25,holder[3],26,holder[3], lty=2)
segments(26,holder[3],26,holder[4], lty=2)



for(i in c(23,25,27)){
  rect(xleft = i, xright = i+1, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("grey", alpha = 0.3))
}  



# Estimate Parameter
library(ergm.count)
b25 <- as.network(y_list_all[[25]], directed=FALSE, matrix.type="a", names.eval="nominations", ignore.eval=FALSE,)
b26 <- as.network(y_list_all[[26]], directed=FALSE, matrix.type="a", names.eval="nominations", ignore.eval=FALSE,)



summary(b25 ~ sum, response="nominations")
summary(b26 ~ sum, response="nominations")



b25.fit <- ergm(b25~sum, response="nominations", reference=~Poisson) 
b26.fit <- ergm(b26~sum, response="nominations", reference=~Poisson) 

summary(b25.fit)
summary(b26.fit)


input_data <- list()
input_data[[1]] <- y_list_all[[25]]; 
input_data[[2]] <- y_list_all[[26]]


# parameter estimation
n <- length(id)
eta <- rep(0, 2)
eta <- partial_stepping_temporal(20, 100, n, input_data, eta, node_attr, edge_attr, pi0=0.2, m=80)
eta <- newton_raphson_temporal(20, 1000, 20*n, input_data, eta, node_attr, edge_attr, pi0=0.2, m=80)

# standard error
se_eta <- SE_temporal(1000, 20*n*n, input_data, eta, node_attr, edge_attr, pi0=0.2, m=80)


eta/ se_eta
write.csv(eta, 'eta.csv', row.names = F)
write.csv(se_eta, "se_eta.csv", row.names = F)


# Support Evidence
counter <- c(0,0)
inc_holder <- dec_holder <- c()

for(i in 1:13){
  for(j in 1:13){
    if(i < j){
      
      if(y_list_all[[26]][i,j] > y_list_all[[25]][i,j]){
        counter[1] <- counter[1] + 1
        inc_holder <- c(inc_holder, y_list_all[[26]][i,j] - y_list_all[[25]][i,j])
      }
      
      if(y_list_all[[26]][i,j] < y_list_all[[25]][i,j]){
        counter[2] <- counter[2] + 1
        dec_holder <- c(dec_holder, y_list_all[[25]][i,j] - y_list_all[[26]][i,j])
      }
      
    }
  }
}


counter
sum(inc_holder)
sum(dec_holder)



