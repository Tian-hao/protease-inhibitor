#Rcode
a1 <- read.table('analysis/potts_fitness.txt',header=T)
png('figures/potts_fitness.png',res=600,height=2800,width=2800)
par(mar=c(4,4,4,4))
plot(a1[,2],a1[,3],xlab='Energy',ylab='fitness',pch=16)
abline(lm(a1[,3]~a1[,2]),lty=2,lwd=2,col='grey')
segments(a1[,2],a1[,3]-a1[,4],a1[,2],a1[,3]+a1[,4],lwd=2)
dev.off()
