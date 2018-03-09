##################################################################
##################################################################
#Get 6 month data and combine with the baseline data
##################################################################
##################################################################
setwd("/Users/Selene/Desktop/6MonthData")

#read raw data
valid=read.csv("/Users/Selene/Desktop/6MonthData/Accel_Only/RfH_6M_WTV_forLoki_20160526.csv",h=T)
data.12=read.csv("/Users/Selene/Desktop/6MonthData/Accel_Only/A2_RfH_2012_AccOnly.csv",h=T)
data.13=read.csv("/Users/Selene/Desktop/6MonthData/Accel_Only/A2_RfH_2013_AccOnly.csv",h=T)
data.14=read.csv("/Users/Selene/Desktop/6MonthData/Accel_Only/A2_RfH_2014_AccOnly.csv",h=T)
data.15=read.csv("/Users/Selene/Desktop/6MonthData/Accel_Only/A2_RfH_2015_AccOnly.csv",h=T)

valid=valid[valid$Valid600==1,]
save(valid,file='valid.rdata')
save(data.12,file='data.12.rdata')
save(data.13,file='data.13.rdata')
save(data.14,file='data.14.rdata')
save(data.15,file='data.15.rdata')

data=rbind(data.12,data.13,data.14,data.15)
save(data,file='data.rdata')
data[[2]]=as.character(data[[2]])
data$date=substr(data[[2]],1,10)
data$time=substr(data[[2]],12,19)

valid[[1]]=as.character(valid[[1]])
valid[[5]]=as.character(valid[[5]])
library(plyr)
list=list()
for(i in 1:length(unique(valid$identifier))){
	x=data[data$identifier==unique(valid$identifier)[i],]
	y=valid[valid$identifier==unique(valid$identifier)[i],]
	l=list()
	for(j in 1:length(unique(y$wearDate))){
		z=x[x$date==unique(y$wearDate)[j],]
		l[[j]]=z
		print(unique(y$wearDate)[j])
	}
	li=ldply(l,data.frame)
	list[[i]]=li

}
fulldat=ldply(list,data.frame)
save(fulldat,file='fulldat.rdata')

fulldat$wearing=as.character(fulldat$wearing)
fulldat$activity=ifelse(fulldat$wearing=='nw',NA,fulldat$activity)
data.f=fulldat[,c('identifier','date','activity','time')]
data.f$identifier=substr(data.f$identifier,4,15)
names(data.f)=c('identifier','dt','activity','time')


#complete times of day (some days don't have the full 1440min record, fill the missing times with activity=NA)
mdat=aggregate(data.f$activity,list(data.f$dt,data.f$identifier),mean,na.rm=TRUE)
names(mdat)=c('dt','identifier','activity')

t1=head(data.f,1440)[,4]
time=rep(t1,nrow(mdat))

l=list()
for(i in 1:length(unique(mdat$identifier))){
	n=sum(mdat$identifier==unique(mdat$identifier)[i])
	l[[i]]=rep(unique(mdat$identifier)[i],n*1440)
}
identifier=unlist(l)

dt=rep(mdat$dt,each=1440)

lis=list()
for(i in 1:length(unique(mdat$identifier))){
	x=data.f[data.f$identifier==unique(mdat$identifier)[i],]
	list=list()
	for(j in 1:length(unique(x$dt))){
		match=match(t1,x[x$dt==unique(x$dt)[j],4])
		match[!is.na(match)]=x[x$dt==unique(x$dt)[j],3]
		list[[j]]=match
	}
	lis[[i]]=unlist(list)
}
activity=unlist(lis)

data.e=data.frame(identifier,dt,activity,time)
save(data.e,file='extended data with activity and date.rdata')


#get proportion missing data for every day
edata=data.e
agg=aggregate(edata$activity,list(edata$dt,edata$identifier),mean,na.rm=TRUE)
names(agg)=c('dt','identifier','activity')
pro=list()
for(i in 1:length(unique(agg$identifier))){
	x=edata[edata$identifier==unique(agg$identifier)[i],]
	s=rep(0,length(unique(x$dt)))
	for(j in 1:length(unique(x$dt))){
		s[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<0,3])/length(x[x$dt==unique(x$dt)[j],3])
	}
	pro[[i]]=s
}
prop=unlist(pro)
agg$prop=prop
ag=na.omit(agg)
save(ag,file='aggregate on mean.rdata')


#get days with >50% data
sub=ag[ag$prop<0.5,]
list=list()
for (i in 1:length(unique(sub$identifier))){
	x=edata[edata$identifier==unique(sub$identifier)[i],]	
	y=sub[sub$identifier==unique(sub$identifier)[i],]
	l=list()
	for(j in 1:length(unique(y$dt))){
		z=x[x$dt==unique(y$dt)[j],]
		l[[j]]=z
	}
	li=ldply(l,data.frame)
	list[[i]]=li
}
newd=ldply(list,data.frame)
save(newd,file='days with more than 50% time.rdata')


#get complete consecutive profiles at the 6 month mark
sub=ag[ag$prop<0.5,]
miss=ifelse(newd$activity=='NA',0,1)
miss[is.na(miss)]=0
newd$miss=miss

list0=list()
for (i in 1:nrow(sub)){
	x=newd[(i-1)*1440+1:1440,]
	r=rle(x$miss)
	list0[[i]]=r
}

s=rep(0,nrow(sub))
for (i in 1:nrow(sub)){
	vec=list0[[i]][1][[1]]
	s[i]=max(vec)
}

whi=which(s<720)
sub=sub[-whi,]
miss=ifelse(edata$activity=='NA',0,1)
miss[is.na(miss)]=0
edata$miss=miss
library(plyr)
list=list()
for (i in 1:length(unique(sub$identifier))){
	x=edata[edata$identifier==unique(sub$identifier)[i],]	
	y=sub[sub$identifier==unique(sub$identifier)[i],]
	l=list()
	for(j in 1:length(unique(y$dt))){
		z=x[x$dt==unique(y$dt)[j],]
		l[[j]]=z
	}
	li=ldply(l,data.frame)
	list[[i]]=li
}
newd=ldply(list,data.frame)
#redo missing patterns + curve registration
list0=list()
for (i in 1:nrow(sub)){
	x=newd[(i-1)*1440+1:1440,]
	r=rle(x$miss)
	list0[[i]]=r
}

s=rep(0,nrow(sub))
for (i in 1:nrow(sub)){
	vec=list0[[i]][1][[1]]
	s[i]=max(vec)
}

list1=list()
for (i in 1:nrow(sub)){
	vec=list0[[i]][1][[1]]
	whic=which(vec==max(vec))
	sum=ifelse(whic==1,0,sum(vec[1:(whic-1)]))
	index=(sum+1):(sum+720)
	x=newd[(i-1)*1440+1:1440,]
	list1[[i]]=x[index,]	
}
gooddf=ldply(list1,data.frame)

save(gooddf,file='complete consecutive profiles.rdata')


#combine baseline data and 6 month data
load("/Users/Selene/Desktop/6MonthData/complete consecutive profiles.rdata")
gooddf$phase=1
gooddf.1=gooddf
gooddf.1$identifier=as.character(gooddf.1$identifier)

load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles.rdata")
gooddf$phase=0
gooddf.0=gooddf
gooddf.0$identifier=substr(gooddf.0$identifier,4,15)

#plot mean curves
y=matrix(gooddf.0$activity,nrow=nrow(gooddf.0)/720,byrow=TRUE)
mu.0=apply(y, 2,mean)
mu.0.sm=ksmooth(1:720,mu.0,kernel='normal',bandwidth=30)
y=matrix(gooddf.1$activity,nrow=nrow(gooddf.1)/720,byrow=TRUE)
mu.1=apply(y, 2,mean)
mu.1.sm=ksmooth(1:720,mu.1,kernel='normal',bandwidth=30)

plot(mu.1.sm,type='l',xlab='time',ylab='mean activity')
lines(mu.0.sm,col='red')
legend("topright",col=c('red','black'),legend=c('baseline','6month'),text.col=c('red','black'),lty=c(1,1))

gooddf.c=rbind(gooddf.0,gooddf.1)
save(gooddf.c,file='combined 2 phase data.rdata')




















