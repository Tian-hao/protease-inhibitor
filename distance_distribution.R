#Rcode
library(RColorBrewer)
mycol <- brewer.pal(4,'Accent')
makeT<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
mycol2 <- makeT(mycol)
a1 <- read.table('analysis/dist_category.txt',header=T)
b1 <- a1[which(a1[,4] == 1),]
b2 <- a1[which(a1[,4] == 2),]
b3 <- a1[which(a1[,4] == 3),]
png('figures/distance_distribution.png',res=600,width=3600,height=3600)
plot(c(0,1),c(0,1),type='n',xlim=c(0,40),ylim=c(0,0.09),xlab='pairwise distance',ylab='density')
polygon(density(b1[,3]),col=mycol2[1],border=mycol[1],lwd=2)
polygon(density(b2[,3]),col=mycol2[2],border=mycol[2],lwd=2)
polygon(density(b3[,3]),col=mycol2[3],border=mycol[3],lwd=2)
legend('topright',legend=c('within DRAM','DRAM and others','among others'),col=mycol,pch=16,bty='n',cex=1.6)
dev.off()
