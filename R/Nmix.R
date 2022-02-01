
Nmix<-
function(y,tag="",seed=0,nsweep=10000,nburnin=0,
	kinit=1,
	qempty=1,qprior=0,qunif=0,qfix=0,qrkpos=0,qrange=1,qkappa=0,qbeta=1,qkreg=0,
      out="",qfull,#qpalloc=0,qpwms=0,qpclass=0,
	alpha=2,beta=0.02,delta=1,eee=0,fff=0,ggg=0.2,hhh=10,unhw=1.0,     
	kappa=1.0,lambda=-1,xi=0.0,sp=1,
	nspace=nsweep%/%1000,nparsamp=7*(nsweep%/%nspace),nsamp=100,nskdel=2,
	nmax=length(y),ncmax=30,ncmax2=10,ncd=7,ngrid=200,kkzz=35,
	idebug=-1,qdebug=0) 
{
if(is.character(y)) {tag<-y; if(exists(y,envir=.GlobalEnv)) y<-get(y,envir=.GlobalEnv) else y<-scan(paste0(y,".dat"))[-1]}
seed<-sdrni(seed) #; cat(seed,'\n')

moves<-'wpah'
if(!qfix) moves<-paste0('s',moves)
if(qempty) moves<-paste0(moves,'b')
imoves<-match(strsplit(moves,NULL)[[1]],strsplit('scawphbd',NULL)[[1]])
# cat(moves,'\n')

aout<-strsplit('dacw','')[[1]]%in%strsplit(out,'')[[1]]
if(any(aout[2:4])) warning("output options a,c,w not yet implemented")
# qpalloc<-aout[2]; qpclass<-aout[3]; qpwms<-aout[4]
qpalloc<-qpclass<-qpwms<-0
if(missing(qfull)||qfull==0) qfull<-any(aout)
#if(qfull) {
ncmax1<-ncmax; ncd1<-ncd; ngrid1<-ngrid; ncmax21<-ncmax2
#} else {ncmax1<-1; ncd1<-1; ngrid1<-1; ncmax21<-1}

z<-.Fortran("nmixsub",
	y=as.single(y),n=as.integer(length(y)),
	post=single(ncmax),den=single(ncd*ngrid),avn=single(ncmax21*ncmax21),
	as.integer(ncmax1),as.integer(ncd1),as.integer(ngrid1),as.integer(ncmax21),
	wtav=single(ncmax2*ncmax2),muav=single(ncmax2*ncmax2),sigav=single(ncmax2*ncmax2),
	acctry=integer(8),enttr=single(nsweep%/%nspace),ktr=integer(nsweep%/%nspace),
	off=integer(nsweep%/%nspace),partr=single(3*nparsamp),as.integer(nparsamp),
	devtr=single(nsweep%/%nspace),
	as.integer(nmax),as.integer(ncmax),as.integer(ncmax2),as.integer(ncd),as.integer(ngrid),as.integer(kkzz),
	nsweep=as.integer(nsweep),nburnin=as.integer(nburnin),as.integer(imoves),as.integer(length(imoves)),
	as.integer(kinit),
	as.integer(qempty),as.integer(qprior),as.integer(qunif),as.integer(qfix),as.integer(qfull),
	as.integer(qrkpos),as.integer(qrange),as.integer(qkappa),as.integer(qbeta),as.integer(qpalloc),
	as.integer(qpwms),as.integer(qpclass),as.integer(qkreg),
	as.single(alpha),as.single(delta),as.single(eee),as.single(fff),as.single(ggg),as.single(hhh),as.single(unhw),
	as.single(kappa),as.single(lambda),as.single(xi),as.integer(sp),
	as.integer(nspace),as.integer(nsamp),as.integer(nskdel),as.integer(idebug),as.single(beta),as.integer(qdebug),
	PACKAGE="Nmix")

yrange<-max(y)-min(y); yd0 = min(y)-0.05*yrange; ydinc = 1.1*yrange/ngrid
if(qfull) {z$den<-matrix(z$den,ncd,ngrid); z$den<-rbind(yd0+ydinc*(1:ngrid),z$den)} else z$den<-NULL

if(!qfull) z$devtr<-NULL
z$traces<-list(k=z$ktr,entropy=z$enttr,deviance=z$devtr)
z$ktr<-NULL; z$enttr<-NULL; z$devtr<-NULL

z$partr<-matrix(z$partr,ncol=3,byrow=TRUE,dimnames=list(NULL,c("wt","mu","sigma")))

zz<-list()
for(t in seq_along(z$off)) 
	{
	tzm<-z$off[t]+z$traces$k[t]
	if(tzm<=nrow(z$partr)) zz[[t]]<-z$partr[z$off[t]+(1:z$traces$k[t]),]
	}
z$partr<-zz
if(length(z$partr)<length(z$off)) warning('parameter sample trace truncated: increase nparsamp')

attr(z$post,'Csingle')<-NULL
z$tag<-tag
z$seed<-seed
z$moves<-moves
z$acctry<-matrix(z$acctry,2,4)
if(qfull) z$avn<-matrix(z$avn,ncmax2,ncmax2) else z$avn<-NULL
dimnames(z$acctry)<-list(c('acc','try'),c('split','merge','birth','death'))

z$wtav<-matrix(z$wtav,ncmax2,ncmax2)
z$muav<-matrix(z$muav,ncmax2,ncmax2)
z$sigav<-matrix(z$sigav,ncmax2,ncmax2)
z$pe<-list()
for(k in 1:ncmax2) if(z$post[k]>0) z$pe[[k]]<-cbind(wt=z$wtav[k,1:k],mu=z$muav[k,1:k],sigma=z$sigav[k,1:k])
z$wtav<-NULL; z$muav<-NULL; z$sigav<-NULL

z[names(z)==""]<-NULL
class(z)<-'nmix'
z
}

