#Rcode
library(RColorBrewer)
mycol <- brewer.pal(4,'Accent')
maket<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
  blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
mycol <- maket(mycol)
mycol[1] <- 1
a1 <- read.table('analysis/fitness_single.txt',header=T)
a1 <- cbind(a1,seq(1,13))
png('figures/fitness_single.png',res=600,height=2400,width=5000)
plot(a1$mean,ylim=c(-3,0.5),type='n',xaxt='n',xlab='',ylab='log10 fitness')
points(jitter(a1[,7]),a1$nd_R1,pch=16,col=mycol[1])
points(jitter(a1[,7]),a1$nd_R2,pch=16,col=mycol[1])
points(jitter(a1[,7]),a1$nd_R3,pch=16,col=mycol[1])
segments(a1[,7],a1$mean-a1$std,a1[,7],a1$mean+a1$std,col=mycol[1],lwd=2)
segments(a1[,7]-0.2,a1$mean-a1$std,a1[,7]+0.2,a1$mean-a1$std,col=mycol[1],lwd=2)
segments(a1[,7]-0.2,a1$mean+a1$std,a1[,7]+0.2,a1$mean+a1$std,col=mycol[1],lwd=2)
axis(side=1,at=seq(1,13),labels=a1$mutation)
dev.off()
