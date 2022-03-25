plplugin<-function(z)
{
opars<-par(no.readonly=TRUE)
on.exit(par(opars))
if(is.null(z$pe)||is.null(z$post)) return()
if(is.null(z$den))
{
if(is.null(z$y)) return()
ngrid<-200
yrange<-max(z$y)-min(z$y); yd0 = min(z$y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
yg<-yd0+ydinc*(1:ngrid)
} else {yg<-z$den[1,]}
layout(1)
par(mar=c(4,4,3,4)+.1)
fyg<-rep(0,length(yg))
for(k in 1:length(z$pe)) if(z$post[k]>0)
{
fyg<-fyg+z$post[k]*apply(z$pe[[k]][,'wt']*dnorm(outer(z$pe[[k]][,'mu'],yg,"-")/z$pe[[k]][,'sigma']),2,sum)
}
lines(yg,fyg,type='l',col='blue',lwd=2)
}
