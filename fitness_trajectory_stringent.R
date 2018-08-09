#Rcode
a1 <- read.table('analysis/stringent_trajectory_fitness.txt',header=T)
png('figures/stringent_trajectory_fitness.png',res=600,height=3000,width=3000)
par(mar=c(4,4,4,4))
plot(c(0,4),c(-4,0),type='n',xlab='Hamming distance to WT',ylab='fitness')
segments(a1[,1],a1[,2],a1[,3],a1[,4],lwd=0.5,col=rgb(0,0,0,alpha=0.5))
dev.off()
