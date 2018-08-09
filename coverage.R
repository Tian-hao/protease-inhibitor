#Rcode
library(RColorBrewer)
mycol <- brewer.pal(8,'Accent')
mycol[4] <- mycol[5]
a1 <- read.table('analysis/coverage.txt',header=T,row.names=1)
a2 <- t(a1)
png('figures/coverage.png',res=600,width=2400,height=2400)
plot(a2[1:11,2]/a2[1:11,1],type='l',lwd=3,col=mycol[1],ylab='coverage',xlab='order of mutations')
points(a2[1:11,3]/a2[1:11,1],type='l',lwd=3,col=mycol[2])
points(a2[1:11,4]/a2[1:11,1],type='l',lwd=3,col=mycol[3])
points(a2[1:11,5]/a2[1:11,1],type='l',lwd=3,col=mycol[4])
legend('topright',legend=c('replicate1','replicate2','replicate3','overlap'),pch=16,col=mycol,bty='n',cex=1)
dev.off()
