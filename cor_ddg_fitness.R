#Rcode
source('pythoncode/fitness_single.R')
mycol <- brewer.pal(6,'Accent')
mycol[4] <- mycol[5]
mycol[5] <- 'black'
mycol <- maket(mycol)
a1 <- read.table('analysis/correlation_ddg_fitness.txt',header=T)
png('figures/correlation_ddg_fitness.png',res=600,width=2400,height=2400)
rho <- round(cor(a1$fitness,a1$ddg,method='spearman'),4)
par(mar=c(4,4,4,4))
plot(a1$fitness,a1$ddg,pch=16,xlab='fitness',ylab='ddg',main=paste('rho=',rho,sep=''),col=mycol[5])
#segments(a1$fitness-a1$fitness_std,a1$ddg,a1$fitness+a1$fitness_std,a1$ddg,lwd=2)
dev.off()

mycol <- brewer.pal(6,'Accent')
mycol[4] <- mycol[5]
png('figures/correlation_ddg_fitness_one2four.png',res=600,width=3600,height=3600)
b1 <- a1[which(a1$mut_number==1),]
b2 <- a1[which(a1$mut_number==2),]
b3 <- a1[which(a1$mut_number==3),]
b4 <- a1[which(a1$mut_number==4),]
rho1 <- round(cor(b1$fitness,b1$ddg,method='spearman'),4)
rho2 <- round(cor(b2$fitness,b2$ddg,method='spearman'),4)
rho3 <- round(cor(b3$fitness,b3$ddg,method='spearman'),4)
rho4 <- round(cor(b4$fitness,b4$ddg,method='spearman'),4)
par(mfrow=c(2,2))
plot(b1$fitness,b1$ddg,pch=16,xlab='fitness',ylab='ddg',main=paste('single mutant\nrho=',rho1,sep=''),col=mycol[1])
segments(b1$fitness-b1$fitness_std,b1$ddg,b1$fitness+b1$fitness_std,b1$ddg,lwd=2,col=mycol[1])

plot(b2$fitness,b2$ddg,pch=16,xlab='fitness',ylab='ddg',main=paste('double mutant\nrho=',rho2,sep=''),col=mycol[2])
segments(b2$fitness-b2$fitness_std,b2$ddg,b2$fitness+b2$fitness_std,b2$ddg,lwd=2,col=mycol[2])

plot(b3$fitness,b3$ddg,pch=16,xlab='fitness',ylab='ddg',main=paste('triple mutant\nrho=',rho3,sep=''),col=mycol[3])
segments(b3$fitness-b3$fitness_std,b3$ddg,b3$fitness+b3$fitness_std,b3$ddg,lwd=2,col=mycol[3])

mycol[4] <- maket(mycol[4])
plot(b4$fitness,b4$ddg,pch=16,xlab='fitness',ylab='ddg',main=paste('quandra mutant\nrho=',rho4,sep=''),col=mycol[4])
segments(b4$fitness-b4$fitness_std,b4$ddg,b4$fitness+b4$fitness_std,b4$ddg,lwd=2,col=mycol[4])
dev.off()
