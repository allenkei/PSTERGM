
obs <- c(2,5,7,4,3,5,9,7,5,5)

plus <- minu <- c()
for(i in 2:length(obs)){
  plus <- c(plus, max(obs[i],obs[i-1]))
  minu <- c(minu, min(obs[i],obs[i-1]))
}

x <- seq(1.5, length(obs)-0.5, by=1)

par(mar=c(5.1, 4.1, 2, 6), xpd=TRUE)
plot(1:length(obs), obs, type='l', xlab='time points', ylab = 'dyad values', 
     lwd = 2, xaxt='n', pch=19, cex = 0.5)


xtick<-1:length(obs)
axis(side=1, at=xtick)
lines(x, plus, col = 'red', type='b', pch=19, lty=2)
lines(x, minu, col = 'blue', type='b', pch=19, lty=2)

lines(c(2.5,3.5), c(7,7), col = 'red', lwd = 3)
lines(c(6.5,7.5), c(9,9), col = 'red', lwd = 3)
lines(c(4.5,5.5), c(3,3), col = 'blue', lwd = 3)


for(i in c(1,3,5,7,9)){
  rect(xleft = i, xright = i+1, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("grey", alpha = 0.3))
}

legend('bottomright', legend=c("Aug", "Obs","Dim"), inset=c(-0.25,0),
       col=c("red", "black","blue"), lty=c(2,1,2), lwd=c(1,2,1), cex=0.8)




