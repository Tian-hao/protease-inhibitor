#Rcode
library(vioplot)
library(RColorBrewer)
mycol <- brewer.pal(4,'Accent')
a1 <- read.table('analysis/fitness_mut_number.txt',header=T)
b1 <- a1[which(a1$MutNumber==1),2]
b2 <- a1[which(a1$MutNumber==2),2]
b3 <- a1[which(a1$MutNumber==3),2]
b4 <- a1[which(a1$MutNumber==4),2]
b5 <- a1[which(a1$MutNumber==5),2]
b6 <- a1[which(a1$MutNumber==6),2]
b7 <- a1[which(a1$MutNumber==7),2]
b8 <- a1[which(a1$MutNumber==8),2]
png('figures/fitness_mut_number.png',res=600,width=3000,height=3000)
vioplot(b1,b2,b3,b4,b5,b6,b7,b8,col=mycol[1])
dev.off()

