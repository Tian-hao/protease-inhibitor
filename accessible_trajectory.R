#Rcode
a1 <- read.table('analysis/accessible_trajectory.txt',header=T)
b2 <- a1[which(a1$hamming_distance==2),]
b3 <- a1[which(a1$hamming_distance==3),]
b4 <- a1[which(a1$hamming_distance==4),]
png('figures/accessible_trajectory.png',res=600,width=6000,height=2400)
par(mfrow=c(1,3))
plot(b2$accessible_frequency,b2$genotype_count,ylim=c(0,max(b2$genotype_count)),xlab='frequency of accessible path',ylab='genotype count',main='HD=2',type='h',lwd=3,lend='square')
plot(b3$accessible_frequency,b3$genotype_count,ylim=c(0,max(b3$genotype_count)),xlab='frequency of accessible path',ylab='genotype count',main='HD=3',type='h',lwd=3,lend='square')
plot(b4$accessible_frequency,b4$genotype_count,ylim=c(0,max(b4$genotype_count)),xlab='frequency of accessible path',ylab='genotype count',main='HD=4',type='h',lwd=3,lend='square')
dev.off()
