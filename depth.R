#Rcode
library(RColorBrewer)
mycol <- brewer.pal(4,'Accent')
mycol <- mycol[1:2]
options(scipen=10)
a1 <- read.table('analysis/depth.txt',header=T,row.names=1)
png('figures/depth.png',res=600,width=3000,height=3600)
mp <- barplot(as.matrix(t(a1[1:6,])),beside=T,main='sequencing depth',xlab='libraries',ylab='readcount',col=mycol,xaxt='n',log='y')
text(mp[1,1:6]+0.5, par('usr')[3], labels = row.names(a1[1:6,]), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
axis(2)
legend('topleft',legend=c('Wildtype','Total depth'),pch=16,col=mycol,bty='n')
dev.off()
