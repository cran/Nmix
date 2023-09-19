plot.nmix<-function(x, which=1:5, separate=FALSE, plugin=FALSE, offset=1, nsamp=50, equi=TRUE, allsort=TRUE, trued=NULL, ...) 
{
opars<-par(no.readonly=TRUE)
on.exit(par(opars))
for(iw in 1:length(which))
{
w<-which[iw]
if(dev.interactive(TRUE)) if(iw>1) if(separate) dev.new() else par(ask=(iw>1)) 
switch(w,
	{
	layout(matrix(c(1, 2), 2, 1), heights = c(2, 1))
	par(mar = c(4, 4, 2, 2))
	pldens(x, plugin, trued)
	np <- length(x$post) + min(0, 3 - match(FALSE, 0 == rev(x$post)))
	barplot(x$post[1:np], names = 1:np)
	title("posterior on k")
	}
,
	pltraces(x$traces,offset)
,
	plsampden(x,nsamp,equi,,trued)
,
	plcoin(x,allsort)
,
	if(exists('pcl',x)) plclass(x)
)
}
invisible(NULL)
}

pldens<-function(x, plugin, trued=NULL)
{
zd<-x$den
if(!is.null(zd))
	{
	zdmax<-max(zd[-1,])
	if(plugin&&!(is.null(x$pe)||is.null(x$post))) 
		{
		yg<-x$den[1,]
		fyg<-rep(0,length(yg))
		for(k in 1:length(x$pe)) if(x$post[k]>0)
		   {
		   w<-x$pe[[k]][,'wt']; m<-x$pe[[k]][,'mu']; s<-x$pe[[k]][,'sigma']
		   fyg<-fyg+x$post[k]*apply(w*dnorm(outer(m,yg,"-")/s)/s,2,sum)
		   }
		zdmax<-max(zdmax,max(fyg))
		}
	if(!is.null(trued)) zdmax<-max(zdmax,max(trued(zd[1,])))
	}
if(!is.null(x$y)) 
	{
	hist(x$y,50,prob=TRUE,main="",xlab="") 
	if(is.null(zd)) title("histogram of data") else 
		{
		if(zdmax>par('usr')[4]) hist(x$y,50,prob=TRUE,main="",xlab="",ylim=c(0,zdmax))
		title("posterior expected density")
		}
	} else {
	if(is.null(zd))
		{
		plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='')
		return(invisible(NULL))
		}
	plot(zd[1,],zd[2,],type='n',ylab='density',xlab='',ylim=c(0,zdmax))
	title("posterior expected density")
	}
mtext(x$tag,1,line=2)
for(j in 2:7) lines(zd[1,],zd[j,])
lines(zd[1,],zd[8,],col='red',lwd=2)
if(!is.null(trued)) lines(zd[1,],trued(zd[1,]),lwd=2,col='blue')
if(!plugin) return()
# add plug-in density estimator
if(is.null(x$pe)||is.null(x$post)) return()
lines(yg,fyg,type='l',col='darkgreen',lwd=2)
}

pltraces<-function(traces,offset=1,ccol=c('black','red','blue','blue'))
{
# only keep simple vectors, not matrices or lists
traces<-traces[names(which(unlist(lapply(traces,function(x){is.vector(x)&!is.list(x)}))))]
if(is.null(traces)||length(traces)==0) return()

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
if(max(zz)-min(zz)<1e-6) {
lab<-signif(mean(zz),6)
u<-rep(0.5,length(zz))
at<-level+0.5
} else {
lab<-signif(pretty(zz,n=3),6)
u<-(zz-min(zz))/(max(zz)-min(zz))
at<-(lab-min(zz))/(max(zz)-min(zz))+level
}
lines(seq_along(u)*ltr/length(u),u+level,col=colr)
axis(side,at=at,labels=lab,col=colr,col.axis=colr)
mtext(names(traces)[itr],side,line=2,at=0.5+level,col=colr)
}
}

plsampden<-function(z,nsamp=50,equi=TRUE,ngrid=200,trued=NULL)
{
if(is.null(z$traces$pars)) return()
layout(1)
opars<-par(mar=c(4,4,3,4)+.1)

yrange<-max(z$y)-min(z$y); yd0 = min(z$y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
yg<-yd0+ydinc*(1:ngrid)

nsweep<-length(z$traces$pars)
if(equi) tw<-1+floor(0.5+((1:nsamp)+runif(1))*nsweep/nsamp)%%nsweep else tw<-sample(nsweep,nsamp)
dd<-matrix(0,ngrid,nsamp)
for(it in seq_along(tw)) for(i in 1:ngrid)
{
t<-tw[it]
wms<-z$traces$pars[[t]]
dd[i,it]<-sum(wms[,1]*dnorm(yg[i],wms[,2],wms[,3]))
}
matplot(yg,dd,type='l',lty=1,ylab='',xlab='',col=c(1,2,3,5,6))
abline(h=0)
mtext("posterior sample of densities",side=3,line=1,cex=1.6)
if(!is.null(trued)) lines(yg,trued(yg),lwd=2,col='blue')
title(xlab=z$tag,line=2)
title(ylab='density',line=2)
}

plcoin<-function(z,allsort=TRUE)
{
if(is.null(z$traces$alloc)) return()
layout(1)
par(mar=c(4,4,3,4)+.1)

n<-z$n
all<-z$traces$alloc
if(allsort) all<-all[order(z$y),]
pa<-matrix(0,n,n)
for(i in 1:n) for(j in 1:n) pa[i,j]<-mean(all[i,]==all[j,])
image(1:n,1:n,pa,asp=1,xlab='',ylab='')
mtext("posterior cluster probabilities",side=3,line=1,cex=1.6)
}

plclass<-function(z,k)
{
opars<-par(no.readonly=TRUE)
on.exit(par(opars))
if(missing(k))
{
#layout(matrix(c(1,3,2,4),2,2))
par(mfrow=c(2,2),mar=c(3,2,3,2)+.1)
o<-order(z$post,decreasing=TRUE)
ks<-sort((o[o<8&o>1])[1:4])
} else {
par(mfrow=c(1,1),mar=c(4,3,3,3)+.1)
ks<-k
}
for(k in ks) if(k>1) 
{
ngrid<-200
yrange<-max(z$y)-min(z$y); yd0 = min(z$y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
yg<-yd0+ydinc*(1:ngrid)
rows<-((k-2)*(k+1))/2+(1:(k-1))
matplot(yg,t(z$pcl[rows,,drop=FALSE]),type='l',lty=1,xlab='',ylab='')
matplot(z$y,t(z$scl[rows,,drop=FALSE]),pch=3,add=TRUE)
mtext('cumulative probability',2,line=2)
mtext(z$tag,1,line=2)
mtext(paste('k =',k),side=3,line=0.5,cex=1.2)
}
}
