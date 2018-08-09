#Rcode
a1 <- read.csv('analysis/coupling_protease_Jpotts.csv',header=F)
png('figures/coupling_CDF.png',res=600,width=2800,height=2800)
par(mar=c(4,4,4,4))
plot(ecdf(a1[,1]),main='',xlab='J',ylab='frequency',verticals=T,pch='',lwd=2)
b1 <- read.csv('analysis/coupling_protease_resistance.csv')
plot(ecdf(b1[,3]),add=T,verticals=T,pch='',col=2,lwd=2)
legend('bottomright',legend=c('all residues','DRAMs'),bty='n',pch=16,col=c(1,2))
dev.off()
