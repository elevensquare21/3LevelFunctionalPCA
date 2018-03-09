##################################################################
##################################################################
#Extend the 2 level FPCA model to 3 level FPCA model
##################################################################
##################################################################
setwd("/Users/Selene/Desktop/6MonthData/2 level fpca")
load("/Users/Selene/Desktop/6MonthData/combined 2 phase data.rdata")
T=720
library(plyr)
l=list()
for(i in 1:n){
	l[[i]]=gooddf.c[gooddf.c$identifier==unique(gooddf.c$identifier)[i],]
}
gooddf=ldply(l,data.frame)
save(gooddf,file='gooddf.rdata')

#visualize mean
n=nrow(gooddf)/720
x=gooddf$activity
y=matrix(x,nrow=n,byrow=TRUE)
t=1:720
mu=apply(y, 2,mean)
plot(mu,type='l')

#perform multilevel functional PCA (3 level to be specific)
resd=matrix(0, nrow=n, ncol=720) 
resd=t( t(y) - mu ) 

index=function(J){
  col.1=gl(J,J-1)
  v=1:J
  col.2=rep(0,J*(J-1))
  for (i in 1:J){
    col.2[((i-1)*(J-1)+1):((i-1)*(J-1)+J-1)]=v[-i]
  }
  cbind(col.1,col.2)
}

index2=function(J1,J2){
	min=min(J1,J2)
	max=max(J1,J2)
	col.1=gl(min,max-1)
	v=1:max
	col.2=rep(0,min*(max-1))
  for (i in 1:min){
    col.2[((i-1)*(max-1)+1):((i-1)*(max-1)+max-1)]=v[-i]
  }
    
  if(J1<=J2){
  	cbind(col.1,col.2)  
   } else{
   	cbind(col.2,col.1)
   }
	
}

#estimate G and Gbd first
gooddf$id=paste(gooddf$identifier,gooddf$phase,sep='_')
M=length(unique(gooddf$id))
J=rep(0,M)
for (i in 1:M){
  J[i]=nrow(gooddf[gooddf$id==unique(gooddf$id)[i],])/720
}
Jsum=J*(J-1)
SUM=sum(J*(J-1))
save(J,file='J.rdata')

mat1 <- matrix(0, SUM, ncol=720)
mat2 <- matrix(0, SUM, ncol=720)
for (m in 1:M){
  if (J[m] > 1){
  for (k in 1:Jsum[m]){
    mat1[ ifelse(m==1,0,sum(Jsum[1:(m-1)]))+k, ]=resd[ifelse(m==1,0,sum(J[1:(m-1)])) + index(J[m])[k,1], ]
     mat2[ ifelse(m==1,0,sum(Jsum[1:(m-1)]))+k, ]=resd[ifelse(m==1,0,sum(J[1:(m-1)])) + index(J[m])[k,2], ]
  }
  }
}

N=720
G <- matrix(0, N, N)
Gbd <- matrix(0, N, N)
for(i in 1:N){ 
  for(j in i:N) {
    G[i,j] <- cov(resd[,i],resd[,j])
    G[j,i] <- G[i,j]
    Gbd[i,j] <- cov(mat1[,i],mat2[,j]) 
    Gbd[j,i] <- Gbd[i,j]
  }
}

save(G,file='G.rdata')
save(Gbd,file='Gbd.rdata')

#estimate Gbp
load("/Users/Selene/Desktop/6MonthData/2 level fpca/gooddf.rdata")
M=length(unique(gooddf$id))
Jmat=matrix(nrow=M,ncol=2)
for (i in 1:M){
  Jmat[i,1]=nrow(gooddf[gooddf$id==unique(gooddf$id)[i] & gooddf$phase==0,])/720
   Jmat[i,2]=nrow(gooddf[gooddf$id==unique(gooddf$id)[i] & gooddf$phase==1,])/720
}
save(Jmat,file='Jmat.rdata')

lis1=list()
lis2=list()
library(plyr)
for(m in 1:M){
	if(Jmat[m,1]>0 & Jmat[m,2]>0 & sum(Jmat[m,])>2){
		ind1=index2(Jmat[m,1],Jmat[m,2])
#		ind2=index2(Jmat[m,2],Jmat[m,1])
		mat1=matrix(0, 2*nrow(ind1),ncol=720)
		mat2=matrix(0, 2*nrow(ind1),ncol=720)
		for(k in 1:  nrow(ind1)){
			mat1[k,]=resd[ifelse(m==1,0,sum(Jmat[1:(m-1),])) + ind1[k,1], ]
			mat2[k,]=resd[ifelse(m==1,0,sum(Jmat[1:(m-1),])) +Jmat[m,1]+ ind1[k,2], ]
		}
		for(k in (nrow(ind1)+1): (2*nrow(ind1))){
			mat1[k,]=resd[ifelse(m==1,0,sum(Jmat[1:(m-1),])) + Jmat[m,1]+ ind1[k-nrow(ind1),2], ]
			mat2[k,]=resd[ifelse(m==1,0,sum(Jmat[1:(m-1),])) + ind1[k-nrow(ind1),1], ]
		}
		lis1[[m]]=mat1
		lis2[[m]]=mat2
	}
}
m1=ldply(lis1,data.frame)
m2=ldply(lis2,data.frame)

N=720
Gbp <- matrix(0, N, N)
for(i in 1:N){ 
  for(j in i:N) {
    Gbp[i,j] <- cov(m1[,i],m2[,j]) 
    Gbp[j,i] <- Gbp[i,j]
  }
}
save(Gbp,file='Gbp.rdata')

#eigen-decomposition of the covariance matrices
N=720
e1 <- eigen(Gbp)
e2 <- eigen(Gbd-Gbp)
e3 <- eigen(G-Gbd)
fpca1.value <- e1$values 
fpca2.value <- e2$values 
fpca3.value <- e3$values 
fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
fpca3.value <- ifelse(fpca3.value>=0, fpca3.value, 0)
percent1 <- (fpca1.value)/sum(fpca1.value)
percent2 <- (fpca2.value)/sum(fpca2.value)
percent3 <- (fpca3.value)/sum(fpca3.value)
K1 <- max( which(cumsum(percent1) < 0.9 | percent1 > 1/N ) + 1)
K2 <- max( which(cumsum(percent2) < 0.9 | percent2 > 1/N ) + 1)
K3 <- max( which(cumsum(percent3) < 0.9 | percent3 > 1/N ) + 1)
rho.s=sum(fpca1.value)/(sum(fpca1.value)+sum(fpca2.value)+sum(fpca3.value))
rho.p=sum(fpca2.value)/(sum(fpca1.value)+sum(fpca2.value)+sum(fpca3.value))
#K1=96 K2=118 K3=312 rho.s=0.095 rho.p=0.078

fpca1.vectors <- e1$vectors[, 1:K1]
fpca2.vectors <- e2$vectors[, 1:K2]
fpca3.vectors <- e3$vectors[, 1:K3]

for(i in 1:K1) {
  v2 <- fpca1.vectors[,i]
  tempsign <- sum(v2)
  fpca1.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
for(i in 1:K2) {
  v2 <- fpca2.vectors[,i]
  tempsign <- sum(v2)
  fpca2.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
for(i in 1:K3) {
  v2 <- fpca3.vectors[,i]
  tempsign <- sum(v2)
  fpca3.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}

save(fpca1.vectors,file='fpca1.vectors.rdata')
save(fpca2.vectors,file='fpca2.vectors.rdata')
save(fpca3.vectors,file='fpca3.vectors.rdata')


#plot first 4 level 1 pca (smoothed)
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca1.vectors[,1]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plus1=mu.sm$y+pc1.sm$y
minus1=mu.sm$y-pc1.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (31.4%)',ylim=c(200,450))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca1.vectors[,2]*1000
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plus2=mu.sm$y+pc2.sm$y
minus2=mu.sm$y-pc2.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (10.3%)',ylim=c(200,450))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca1.vectors[,3]*1000
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plus3=mu.sm$y+pc3.sm$y
minus3=mu.sm$y-pc3.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (7.4%)',ylim=c(200,450))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca1.vectors[,4]*1000
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plus4=mu.sm$y+pc4.sm$y
minus4=mu.sm$y-pc4.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (4.0%)',ylim=c(200,450))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')

#plot first 4 level 2 pca (smoothed)
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca2.vectors[,1]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plus1=mu.sm$y+pc1.sm$y
minus1=mu.sm$y-pc1.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (11.9%)',ylim=c(200,450))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca2.vectors[,2]*1000
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plus2=mu.sm$y+pc2.sm$y
minus2=mu.sm$y-pc2.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (7.6%)',ylim=c(200,450))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca2.vectors[,3]*1000
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plus3=mu.sm$y+pc3.sm$y
minus3=mu.sm$y-pc3.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (5.9%)',ylim=c(200,450))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca2.vectors[,4]*1000
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plus4=mu.sm$y+pc4.sm$y
minus4=mu.sm$y-pc4.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (4.0%)',ylim=c(200,450))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')


#calculate principal component scores using the projection method introduced in Di's paper
cross.integral.C <- matrix(0, K1, K2)
for(i in 1:K1)
  for(j in 1:K2) 
    cross.integral.C[i,j] <- sum(fpca1.vectors[,i]* fpca2.vectors[,j]) 

cross.integral.D <- matrix(0, K1, K3)
for(i in 1:K1)
  for(j in 1:K3) 
    cross.integral.D[i,j] <- sum(fpca1.vectors[,i]* fpca3.vectors[,j]) 

cross.integral.E <- matrix(0, K2, K3)
for(i in 1:K2)
  for(j in 1:K3) 
    cross.integral.E[i,j] <- sum(fpca2.vectors[,i]* fpca3.vectors[,j]) 

dim(cross.integral.C)
dim(cross.integral.D)
dim(cross.integral.E)

n=nrow(gooddf)/720
int1 <- matrix(0, n, K1)
int2 <- matrix(0, n, K2)
int3 <- matrix(0, n, K3)
for(i in 1:n)   {
  for(j in 1:K1)  {
    int1[ i ,j] <- sum( resd[i,] * fpca1.vectors[,j] ) 
  }
  for(j in 1:K2) {
    int2[ i ,j] <- sum( resd[i,] * fpca2.vectors[,j] )    
  }
  for(j in 1:K3) {
    int3[ i ,j] <- sum( resd[i,] * fpca3.vectors[,j] )    
  }
}

dim(int1)
dim(int2)
dim(int3)

s1 <- matrix(0, n, K1)
s2 <- matrix(0, n, K2)
s3 <- matrix(0, n, K3)
library(MASS)

IminusEE=ginv(diag(rep(1,K2))-cross.integral.E %*% t(cross.integral.E))
CtminusED=t(cross.integral.C)-cross.integral.E %*% t(cross.integral.D)
CminusDE=cross.integral.C-cross.integral.D %*% t(cross.integral.E)
design.xi <- ginv(diag(rep(1,K1))-cross.integral.D %*% t(cross.integral.D)-CminusDE %*% IminusEE %*% CtminusED)

M=length(unique(gooddf$identifier))
J=rep(0,M)
for(i in 1:M){
	J[i]=sum(Jmat[i,])
}
for(m in 1:M){
	resid<- rep(0, K1)
	for(j in 1:J[m]) {
    index <-  ifelse(m==1,0,sum(J[1:(m-1)])) + j
    resid <- resid + ( int1[index,] - cross.integral.D %*% int3[index,]-CminusDE %*% IminusEE %*% (int2[index,]-cross.integral.E %*% int3[index,]) )/J[m]
    }
    index.m <- ( ifelse(m==1,0,sum(J[1:(m-1)])) + 1 ) : (sum(J[1:m]))
    xi.temp <- design.xi %*% resid
    s1[index.m,] <- matrix(rep(xi.temp, each=J[m]), nrow=J[m])


}

load("/Users/Selene/Desktop/6MonthData/2 level fpca/J.rdata")
Jp=J
M=length(Jp)
for(m in 1:M){
	resid<- rep(0, K2)
	for(j in 1:Jp[m]) {
    index <-  ifelse(m==1,0,sum(Jp[1:(m-1)])) + j
    b1=s1[index,]
    resid <- resid + (int2[index,]-cross.integral.E %*% int3[index,]-CtminusED %*% b1)/Jp[m]
    }
    index.m <- ( ifelse(m==1,0,sum(Jp[1:(m-1)])) + 1 ) : (sum(Jp[1:m]))
    xi.temp <- IminusEE %*% resid
    s2[index.m,] <- matrix(rep(xi.temp, each=Jp[m]), nrow=Jp[m])	

	for(j in 1:Jp[m]) {
    index <-  ifelse(m==1,0,sum(Jp[1:(m-1)])) + j
    b1=s1[index,]
    b2=xi.temp
    s3[index,] <- int3[index,]-t(cross.integral.D) %*% b1 -t(cross.integral.E) %*% b2
    		
	}
			
}

coef.1=s1
save(coef.1,file='coef.1.rdata')
coef.2=s2
save(coef.2,file='coef.2.rdata')
coef.3=s3
save(coef.3,file='coef.3.rdata')

















