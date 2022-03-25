readf2cio<-function(fn,imax=Inf,verbose=FALSE)
{
# binary files holding lists, matrices or vectors of numeric data
# writable from Fortran via f2cio interface, readable in R using readBin
# file structure supported: binary file, with integer(4), real(4) or double(8) data
# record 1: list: 0 0
#           matrix or vector: nc mode (mode = 1, 2 or 3 for integer(4), real(4) or double(8))
# succeeding records, one per component of list or row of matrix: 
#           list: number of items, mode as integers, followed by data for this component
#           (note that modes can differ between but not within components)
#           matrix or vector: data for this row
#           matrices of one column are delivered as vectors
#
con<-file(fn,'rb')
na<-readBin(con,'integer',2,size=4)
if(na[1]==0) {
# list case
   if(verbose) cat('reading a list from', fn,'\n')
   res<-list()
   i<-1
   repeat{
	na<-readBin(con,'integer',2,size=4)
	if(length(na)==0) break
	nc<-na[1]
	mode<-c('integer','numeric','double')[na[2]]
	size<-c(4,4,8)[na[2]]
      if(verbose) cat(paste0('component ',i,' is of length ',nc,' values, mode ',mode,' and size ',size,' bytes'),'\n')
	res[[i]]<-readBin(con,mode,nc,size=size)
	i<-i+1
      if(i>imax) break
   }
} else {
# matrix/vector case
   nc<-na[1]
   mode<-c('integer','numeric','double')[na[2]]
   size<-c(4,4,8)[na[2]]
   res<-matrix(0,0,nc)
   repeat{
	row<-readBin(con,mode,nc,size=size)
	if(length(row)==0) break
	res<-rbind(res,row)
   }
   if(nc==1) res<-as.vector(res) else attr(res,'dimnames')<-NULL
}
close(con)
res
}
