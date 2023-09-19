Nmix<-
function(y,tag="",seed=0,nsweep=10000,nburnin=0,
	kinit=1,qempty=1,qprior=0,qunif=0,qfix=0,qrkpos=0,qrange=1,qkappa=0,qbeta=1,
	alpha=2,beta=0.02,delta=1,eee=0,fff=0,ggg=0.2,hhh=10,unhw=1.0,     
	kappa=1.0,lambda=-1,xi=0.0,sp=1,
      out="Dkdep",nspace=nsweep%/%1000,
	nmax=length(y),ncmax=30,ncmax2=10,ncd=7,ngrid=200,k1k2=c(2,8),
	idebug=-1,qdebug=0) 
{
if(is.character(y)) {tag<-y; if(exists(y,envir=.GlobalEnv)) y<-get(y,envir=.GlobalEnv) else y<-scan(paste0(y,".dat"))[-1]}
n<-as.integer(length(y))
seed<-sdrni(seed) #; cat(seed,'\n')

moves<-'wpah'
if(!qfix) moves<-paste0('s',moves)
if((!qfix)&qempty) moves<-paste0(moves,'b')
imoves<-match(strsplit(moves,NULL)[[1]],strsplit('scawphbd',NULL)[[1]])
# cat('moves',moves,imoves,'\n')

if(out=='*') out<-'DCApkdea'
sout<-unlist(strsplit('DCA',''))%in%unlist(strsplit(out,''))
tout<-unlist(strsplit('pkdea',''))%in%unlist(strsplit(out,''))

k1<-k1k2[1]; k2<-k1k2[2]; kkzz<-((k2-k1+1)*(k1+k2))/2

tfnames<-c('partr','ktr','devtr','enttr','ztr')
stem<-paste0(tempdir(),'/',substring(basename(tempfile()),5,16))
tf<-paste0(stem,'.',tfnames)
names(tf)<-tfnames
zfn<-encode(tf); ifn<-zfn[[1]]; lfn<-zfn[[2]]; nfn<-ncol(ifn)

z<-.Fortran("nmixsub",
	y=as.single(y),n=as.integer(n),
	post=single(ncmax),den=single(ncd*ngrid),avn=single(ncmax2*ncmax2),
	wtav=single(ncmax2*ncmax2),muav=single(ncmax2*ncmax2),sigav=single(ncmax2*ncmax2),
	acctry=integer(8),
      as.integer(ifn),as.integer(lfn),as.integer(nfn),
	pcl=single(kkzz*ngrid),scl=single(kkzz*nmax),
	as.integer(nmax),as.integer(ncmax),as.integer(ncmax2),as.integer(ncd),as.integer(ngrid),
	as.integer(k1),as.integer(k2),as.integer(kkzz),
	nsweep=as.integer(nsweep),nburnin=as.integer(nburnin),as.integer(imoves),as.integer(length(imoves)),
	as.integer(kinit),as.integer(nspace),as.integer(idebug),
	as.integer(c(qempty,qprior,qunif,qfix,qrkpos,qrange,qkappa,qbeta,sout,tout,qdebug)),
	as.single(c(alpha,beta,delta,eee,fff,ggg,hhh,unhw,kappa,lambda,xi,sp)),
	iflag=integer(1),
	PACKAGE="Nmix")

# cat('z$post',z$post,'\n')

# traces
partr<-ktr<-devtr<-enttr<-ztr<-NULL
if(tout[1]) {
	partr<-readf2cio(tf['partr'])
	for(i in seq_along(partr)) partr[[i]]<-
		matrix(partr[[i]],ncol=3,byrow=TRUE,dimnames=list(NULL,c("wt","mu","sigma")))
}
if(tout[2]) ktr<-readf2cio(tf['ktr'])
if(tout[3]) devtr<-readf2cio(tf['devtr'])
if(tout[4]) enttr<-readf2cio(tf['enttr'])
if(tout[5]) ztr<-t(readf2cio(tf['ztr']))
z$traces<-list(pars=partr,k=ktr,deviance=devtr,entropy=enttr,alloc=ztr)

# summaries
yrange<-max(y)-min(y); yd0 = min(y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
if(sout[1]) {z$den<-matrix(z$den,ncd,ngrid); z$den<-rbind(yd0+ydinc*(1:ngrid),z$den)} else z$den<-NULL
if(sout[2]) {z$pcl<-matrix(z$pcl,kkzz,ngrid); z$scl<-matrix(z$scl,kkzz,nmax)} else z$pcl<-z$scl<-NULL
if(sout[3]) z$avn<-matrix(z$avn,ncmax2,ncmax2) else z$avn<-NULL

# posterior expectations of parameters
z$wtav<-matrix(z$wtav,ncmax2,ncmax2)
z$muav<-matrix(z$muav,ncmax2,ncmax2)
z$sigav<-matrix(z$sigav,ncmax2,ncmax2)
z$pe<-list()
for(k in 1:ncmax2) if(z$post[k]>0) z$pe[[k]]<-cbind(wt=z$wtav[k,1:k],mu=z$muav[k,1:k],sigma=z$sigav[k,1:k])
z$wtav<-NULL; z$muav<-NULL; z$sigav<-NULL

# miscellaneous
z$stem<-stem
z$call<-match.call()
attr(z$post,'Csingle')<-NULL
z$tag<-tag
z$seed<-seed
z$moves<-moves
z$acctry<-matrix(z$acctry,2,4)
dimnames(z$acctry)<-list(c('acc','try'),c('split','merge','birth','death'))

z[names(z)==""]<-NULL
class(z)<-'nmix'
z
}

