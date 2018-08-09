#Rcode
a1 <- read.table('analysis/validation_fitness_correlation.txt',header=T)
a1 <- a1[which(a1[,4] != 'NA'),]
png('figures/validation_fitness.png',width=4000,height=1800,res=600)
par(mfrow=c(1,3))
plot(a1$fitness,a1$validation,xlim=c(-2.5,0.5),ylim=c(-2.5,0.5),pch=16,xlab='screening',ylab='validation',main=paste('rho=',round(cor(a1$fitness,a1$validation,method='spearman'),4),sep=''))
segments(a1$fitness,a1$validation-a1$validation_std,a1$fitness,a1$validation+a1$validation_std,lwd=2)
segments(a1$fitness-a1$fitness_std,a1$validation,a1$fitness+a1$fitness_std,a1$validation,lwd=2)
abline(lm(a1$validation ~ a1$fitness),lwd=2,lty=2,col='grey')

plot(a1$fitness,a1$mixed_validation,xlim=c(-2.5,0.5),ylim=c(-2.5,0.5),pch=16,xlab='screening',ylab='mixed_validation',main=paste('rho=',round(cor(a1$fitness,a1$mixed_validation,method='spearman'),4),sep=''))
segments(a1$fitness,a1$mixed_validation-a1$mixed_validation_std,a1$fitness,a1$mixed_validation+a1$mixed_validation_std,lwd=2)
segments(a1$fitness-a1$fitness_std,a1$mixed_validation,a1$fitness+a1$fitness_std,a1$mixed_validation,lwd=2)
abline(lm(a1$mixed_validation ~ a1$fitness),lwd=2,lty=2,col='grey')

plot(a1$validation,a1$mixed_validation,xlim=c(-2.5,0.5),ylim=c(-2.5,0.5),pch=16,xlab='validation',ylab='mixed_validation',main=paste('rho=',round(cor(a1$validation,a1$mixed_validation,method='spearman'),4),sep=''))
segments(a1$validation,a1$mixed_validation-a1$mixed_validation_std,a1$validation,a1$mixed_validation+a1$mixed_validation_std,lwd=2)
segments(a1$validation-a1$validation_std,a1$mixed_validation,a1$validation+a1$validation_std,a1$mixed_validation,lwd=2)
abline(lm(a1$mixed_validation ~ a1$validation),lwd=2,lty=2,col='grey')
dev.off()
