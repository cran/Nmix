reimport<-function(z)
{
# reads temporary binary files created in Nmix back into R

#q<-system(paste0("ls ",gsub('\\\\','/',z$stem),"*"),intern=TRUE)
# q<-paste0(tempdir(),'/',list.files(tempdir(),basename(z$stem)))
# or perhaps better dirname(z$stem) instead of tempdir():
q<-paste0(dirname(z$stem),'/',list.files(dirname(z$stem),basename(z$stem)))
res<-list()
for(i in seq_along(q)) res[[i]]<-readf2cio(q[i])
ss<-strsplit(q,".",fixed=TRUE)
exts<-unlist(lapply(ss,function(x){x[length(x)]}))
names(res)<-exts
res
}

# to match what is done in Nmix, would need to change filename extensions to k, pars etc
# pars component is not re-shaped as 3-col matrix, and columns not named 