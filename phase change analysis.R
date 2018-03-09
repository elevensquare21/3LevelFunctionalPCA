##################################################################
##################################################################
#Regression analysis on phase change 
##################################################################
##################################################################
load("/Users/Selene/Desktop/6MonthData/2 level fpca/coef.1.rdata")
load("/Users/Selene/Desktop/6MonthData/2 level fpca/coef.2.rdata")
load("/Users/Selene/Desktop/6MonthData/2 level fpca/coef.3.rdata")
load("/Users/Selene/Desktop/6MonthData/2 level fpca/Jmat.rdata")
load("/Users/Selene/Desktop/6MonthData/2 level fpca/gooddf.rdata")

#get information on which intervention arm each subject is assigned to
nv=read.csv('/Users/Selene/Desktop/ReachForHealth Data/Reach for Health.csv',h=T)
nv=nv[,c(1,16)]
nv=nv[order(nv$ID),]
nv2=read.csv('/Users/Selene/Desktop/ReachForHealth Data/RfH_ScreenCode2ID.csv',h=T)
nv2=nv2[order(nv2$PIDCC),]

newd=data.frame(id=nv2$PtID,arm=nv$ArmAssigned)

match=match(unique(gooddf$identifier),newd[[1]])
armdf=data.frame(id=unique(gooddf$identifier), arm=newd[match,2])

#regression of the difference in each phase coefficient between baseline and 6 month on intervention arm
#save all summary results in a list 
phase.list=list()
for(i in 1:ncol(phase.coef)){
	Jvec=as.vector(t(Jmat))
	Jvec[Jvec!=0]=phase.coef[,i]
	mat=matrix(Jvec,ncol=2,byrow=TRUE)
	mat[mat==0]=NA
	regdf=data.frame(diff=mat[,2]-mat[,1],arm=armdf[,2])
	regdf$arm=as.factor(regdf$arm)
	regdf=na.omit(regdf)
	phase.list[[i]]=summary(lm(diff~arm,data=regdf))

}
setwd('/Users/Selene/Desktop')
save(phase.list,file='phase.list.rdata')

#get the ones that result in significant difference
pvec=rep(NA,ncol(phase.coef))
for(i in 1:ncol(phase.coef)){
	f=phase.list[[i]]$fstatistic
	pvec[i]=pf(f[1],f[2],f[3],lower.tail=F)
}
which(pvec<0.05)

#plot the first four significant ones
load("/Users/Selene/Desktop/6MonthData/2 level fpca/fpca2.vectors.rdata")
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca2.vectors[,which(pvec<0.05)[1]]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plus1=mu.sm$y+pc1.sm$y
minus1=mu.sm$y-pc1.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 6',ylim=c(200,450))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca2.vectors[,which(pvec<0.05)[2]]*3000
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plus2=mu.sm$y+pc2.sm$y
minus2=mu.sm$y-pc2.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 29',ylim=c(200,450))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca2.vectors[,which(pvec<0.05)[3]]*3000
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plus3=mu.sm$y+pc3.sm$y
minus3=mu.sm$y-pc3.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 32',ylim=c(200,450))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca2.vectors[,which(pvec<0.05)[4]]*5000
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plus4=mu.sm$y+pc4.sm$y
minus4=mu.sm$y-pc4.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 61',ylim=c(200,450))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')

#plot mean curves
load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles.rdata")
gooddf$identifier=as.character(gooddf$identifier)
gooddf$identifier=substr(gooddf$identifier,4,15)
n=nrow(gooddf)/720
for(i in 1:n){
	ID=gooddf$identifier[(i-1)*720+1]
	ma=match(ID,armdf$id)
	gooddf$Arm[((i-1)*720+1):(i*720)]=armdf[ma,2]
}
dt11=gooddf[gooddf$Arm==11,]
dt12=gooddf[gooddf$Arm==12,]
dt13=gooddf[gooddf$Arm==13,]
dt14=gooddf[gooddf$Arm==14,]

x.11=dt11$activity
y.11=matrix(x.11,nrow=nrow(dt11)/720,byrow=TRUE)
mu.11=apply(y.11, 2,mean)
sm.11=ksmooth(1:720,mu.11,kernel='normal',bandwidth=30)

x.12=dt12$activity
y.12=matrix(x.12,nrow=nrow(dt12)/720,byrow=TRUE)
mu.12=apply(y.12, 2,mean)
sm.12=ksmooth(1:720,mu.12,kernel='normal',bandwidth=30)

x.13=dt13$activity
y.13=matrix(x.13,nrow=nrow(dt13)/720,byrow=TRUE)
mu.13=apply(y.13, 2,mean)
sm.13=ksmooth(1:720,mu.13,kernel='normal',bandwidth=30)

x.14=dt14$activity
y.14=matrix(x.14,nrow=nrow(dt14)/720,byrow=TRUE)
mu.14=apply(y.14, 2,mean)
sm.14=ksmooth(1:720,mu.14,kernel='normal',bandwidth=30)

plot(sm.11,type='l',xlab='time',ylab='activity counts',lwd=3)
lines(sm.12,lwd=3,col='red')
lines(sm.13,lwd=3,col='blue')
lines(sm.14,lwd=3,col='green')
legend("topright",col=c('black','red','blue','green'),legend=c('11','12','13','14'),text.col=c('black','red','blue','green'),lty=c(1,1,1,1))


















