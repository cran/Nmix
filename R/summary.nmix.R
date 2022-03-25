summary.nmix<-function(object, ...) 
{
cat('Call: '); print(object$call)
print(object)
print(object$acctry)
nzt<-names(object$traces)[!sapply(object$traces,is.null)]
if(length(nzt)==0) nzt<-'none'
cat('Posterior traces available:',nzt,'\n')
invisible(object)
}
