plot.nmix<-function(x, which=1, ...) 
{
if(1%in%which)
{
	oldpar<-par(no.readonly=TRUE)
	on.exit(par(oldpar))
	layout(matrix(c(1, 2), 2, 1), heights = c(2, 1))
	par(mar = c(4, 4, 2, 2),ask=TRUE)
	pldens(x, add = TRUE)
	np <- length(x$post) + min(0, 3 - match(FALSE, 0 == rev(x$post)))
	barplot(x$post[1:np], names = 1:np)
	title("posterior on k")
}
if(2%in%which) pltraces(x$traces)
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

