plot.nmix<-function(x, which=1, offset=1, nsamp=50, equi=TRUE, ...) 
{
opars<-par(no.readonly=TRUE)
on.exit(par(opars))
for(iw in 1:length(which))
{
w<-which[iw]
if(iw>1) par(ask=(iw>1)) 
switch(w,
	{
	layout(matrix(c(1, 2), 2, 1), heights = c(2, 1))
	par(mar = c(4, 4, 2, 2))
	pldens(x, add = TRUE)
	np <- length(x$post) + min(0, 3 - match(FALSE, 0 == rev(x$post)))
	barplot(x$post[1:np], names = 1:np)
	title("posterior on k")
	}
,
	pltraces(x$traces,offset)
,
	plsampden(x,nsamp,equi)
)
}
invisible(NULL)
}

pldens<-function(x='fort.13',add=FALSE)
{
if(is.character(x)) 
	{
	zd<-scan(x) 
	zd<-matrix(zd,8,200,byrow=T)
	} else {
	if(is.list(x)) zd<-x$den else zd<-x
	}
if(add&!is.null(x$y)) 
	{
	hist(x$y,50,prob=TRUE,main="",xlab="") 
	if(is.null(zd)) title("histogram of data") else title("posterior expected density")
	} else {
	if(is.null(zd))
		{
		plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='')
		return(invisible(NULL))
		}
	plot(zd[1,],zd[2,],type='n',ylab='density',xlab='',ylim=c(0,max(zd[2:8,])))
	title("posterior expected density")
	}
mtext(x$tag,1,line=2)
for(j in 2:7) lines(zd[1,],zd[j,])
lines(zd[1,],zd[8,],col='red',lwd=2)
}

pltraces<-function(traces,offset=1,ccol=c('black','red','blue','darkgreen'))
{
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

plsampden<-function(z,nsamp=50,equi=TRUE,ngrid=200)
{
layout(1)
par(mar=c(4,4,3,4)+.1)

yrange<-max(z$y)-min(z$y); yd0 = min(z$y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
yg<-yd0+ydinc*(1:ngrid)

nsweep<-length(z$partr)
if(equi) tw<-1+floor(0.5+((1:nsamp)+runif(1))*nsweep/nsamp)%%nsweep else tw<-sample(nsweep,nsamp)
dd<-matrix(0,ngrid,nsamp)
for(it in seq_along(tw)) for(i in 1:ngrid)
{
t<-tw[it]
dd[i,it]<-sum((dnorm((yg[i]-z$partr[[t]][,2])/z$partr[[t]][,3]))*z$partr[[t]][,1])
}
matplot(yg,dd,type='l',lty=1,ylab='',xlab='')
mtext("Posterior sample",side=3,line=1,cex=1.6)
title(xlab=z$tag,line=2)
title(ylab='density',line=2)
}

