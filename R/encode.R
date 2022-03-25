encode<-function(fn)
{
nr<-length(fn)
ifn = matrix(as.integer(0), nrow = 200, ncol = nr)
  lfn = rep(as.integer(0), length = nr)
  for (i in 1:nr) {
    nchars <- as.integer(nchar(as.character(fn[i])))       
    ifn[1:nchars,i] <- as.integer(unlist(sapply(as.character(fn[i]),charToRaw),use.names=FALSE))
    lfn[i] = nchars
  }
ifn<-ifn[1:max(lfn),,drop=FALSE]
list(ifn,lfn)
}
