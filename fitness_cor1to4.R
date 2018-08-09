#Rcode
library(RColorBrewer)
mycol <- brewer.pal(8,'Accent')
mycol[4] <- mycol[5]
a1 <- read.table('analysis/fitness_one2four.txt',header=T)

drawfigures <- function(mutnum){
if (mutnum <= 3){
  a2 <- a1[which(a1$mut_number==mutnum),]
}
else{
  a2 <- a1[which(a1$mut_number>=mutnum),]
}
a2$mut_number <- mutnum
png(paste('figures/fitness_cor12_',mutnum,'mut.png',sep=''),res=600,height=2400,width=2400)
par(mar=c(4,4,4,4))
plot(log10(a2$nd_R1),log10(a2$nd_R2),col=mycol[a2$mut_number],pch=16,xlab='',ylab='',xlim=c(-4,1),ylim=c(-4,1))
title(ylab='replicate2',line=2,cex.lab=1.5)
title(xlab='replicate1',line=2,cex.lab=1.5)
legend('bottomright',legend=paste('rho=',round(cor(log10(a2$nd_R1),log10(a2$nd_R2),method='spearman'),2),sep=''),bty='n',cex=1.5)
abline(0,1,col='grey',lwd=2,lty=2)
dev.off()

png(paste('figures/fitness_cor13_',mutnum,'mut.png',sep=''),res=600,height=2400,width=2400)
par(mar=c(4,4,4,4))
plot(log10(a2$nd_R1),log10(a2$nd_R3),col=mycol[a2$mut_number],pch=16,xlab='',ylab='',xlim=c(-4,1),ylim=c(-4,1))
title(ylab='replicate3',line=2,cex.lab=1.5)
title(xlab='replicate1',line=2,cex.lab=1.5)
legend('bottomright',legend=paste('rho=',round(cor(log10(a2$nd_R1),log10(a2$nd_R3),method='spearman'),2),sep=''),bty='n',cex=1.5)
abline(0,1,col='grey',lwd=2,lty=2)
dev.off()

png(paste('figures/fitness_cor23_',mutnum,'mut.png',sep=''),res=600,height=2400,width=2400)
par(mar=c(4,4,4,4))
plot(log10(a2$nd_R2),log10(a2$nd_R3),col=mycol[a2$mut_number],pch=16,xlab='',ylab='',xlim=c(-4,1),ylim=c(-4,1))
title(ylab='replicate3',line=2,cex.lab=1.5)
title(xlab='replicate2',line=2,cex.lab=1.5)
legend('bottomright',legend=paste('rho=',round(cor(log10(a2$nd_R2),log10(a2$nd_R3),method='spearman'),2),sep=''),bty='n',cex=1.5)
abline(0,1,col='grey',lwd=2,lty=2)
dev.off()
}

drawfigures(1)
drawfigures(2)
drawfigures(3)
drawfigures(4)
