pltraces<-function(traces,offset=1,ccol=c('black','red','blue','darkgreen'))
{
oldpar<-par(no.readonly=TRUE)
on.exit(par(oldpar))
layout(1)
par(mar=c(4,4,3,4)+.1)

ntr<-length(traces)
ltr<-max(sapply(traces,length))
plot(c(0,ltr),c(0,1+offset*(ntr-1)),type='n',yaxt='n',xlab='',ylab='')
title(xlab='thinned sample',line=2)
for(itr in 1:ntr)
{
colr<-ccol[1+(itr-1)%%length(ccol)]
level<-offset*(ntr-itr)
side<-4-2*(itr%%2)
zz<-traces[[itr]]
lab<-pretty(zz,n=3)
u<-(zz-min(zz))/(max(zz)-min(zz))
lines(seq_along(u)*ltr/length(u),u+level,col=colr)
axis(side,at=(lab-min(zz))/(max(zz)-min(zz))+level,labels=lab,col=colr,col.axis=colr)
mtext(names(traces)[itr],side,line=2,at=0.5+level,col=colr)
}

}
