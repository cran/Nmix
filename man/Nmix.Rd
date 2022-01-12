\name{Nmix}
\alias{Nmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Inference on Univariate Normal Mixtures
}
\description{
Wrapper for Nmix Fortran program
}
\usage{
Nmix(y,tag="",seed=0,nsweep=10000,nburnin=0,kinit=1,
	qempty=1,qprior=0,qunif=0,qfix=0,qrkpos=0,qrange=1,qkappa=0,
	qbeta=1,qkreg=0,out="",qfull,alpha=2,beta=0.02,delta=1,eee=0,
	fff=0,ggg=0.2,hhh=10,unhw=1.0,kappa=1.0,lambda=-1,xi=0.0,sp=1,
	nspace=nsweep\%/\%1000,nsamp=100,nskdel=2,nmax=length(y),ncmax=30,
	ncmax2=10,ncd=7,ngrid=200,kkzz=35,idebug=-1,qdebug=0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
either (i) a numerical data vector, (ii) a character scalar naming a numerical data vector in the global environment 
or (iii) a character scalar naming a dataset to be found in file y.dat
}
  \item{tag}{
name for the dataset, in the case that \code{y} is a numerical vector
}
  \item{seed}{
positive integer to set random number seed for a reproducible run, or 0 to initialise this process; output value
can be used to replicate run subsequently 
}
  \item{nsweep}{
number of sweeps
}
  \item{nburnin}{
length of burn in
}
  \item{nspace}{
spacing between samples recorded in time-series traces (\code{ktr}, etc)
}
  \item{out}{
character string to specify optional output: string containing letters 'a','c','d',or 'w' correspond to Fortran program options -pa, -pc, -pd, or -pw. Only 'd' currently implemented
}
  \item{qfull}{
(deprecated) integer or logical scalar: 1 or TRUE is equivalent to \code{out} containing \code{"d"})
}
  \item{kinit}{
integer, initial number of components
}
  \item{qempty}{
integer, 1 or 0 according to whether the empty-component birth/death moves should be used
}
  \item{qprior}{
integer, 1 or 0 according to whether the prior should be simulated instead of the posterior
}
  \item{qunif}{
integer, 1 or 0 according to whether the uniform proposals should be used for the component means instead of gaussian ones
}
  \item{qfix}{
integer, 1 or 0 according to whether the number of components should be held fixed (at the value of \code{kinit})
}
  \item{qrkpos}{
integer, 1 or 0 according to whether the the number of non-empty components should be reported throughout
}
  \item{qrange}{
integer, 1 or 0 according to whether range-based parameter priors should be used
}
  \item{qkappa}{
integer, 1 or 0 according to whether \code{kappa} should be updated
}
  \item{qbeta}{
integer, 1 or 0 according to whether \code{beta} should be updated
}
  \item{qkreg}{
integer, 1 or 0 according to whether
}
  \item{alpha}{
numeric, set value of parameter alpha
}
  \item{beta}{
numeric, set value of parameter beta
}
  \item{delta}{
numeric, set value of parameter delta
}
  \item{eee}{
numeric, set value of parameter e
}
  \item{fff}{
numeric, set value of parameter f
}
  \item{ggg}{
numeric, set value of parameter g
}
  \item{hhh}{
numeric, set value of parameter h
}
  \item{unhw}{
numeric, set value of half-width for uniform proposals
}
  \item{kappa}{
numeric, set value of parameter kappa
}
  \item{lambda}{
numeric, set value of parameter lambda
}
  \item{xi}{
numeric, set value of parameter xi
}
  \item{sp}{
numeric, set value of parameter s
}
  \item{nsamp}{
integer, set thinning interval for some outputs
}
  \item{nskdel}{
integer, set 
}
  \item{nmax}{
integer, set upper bound for \code{n}
}
  \item{ncmax}{
integer, set upper bound for \code{k}
}
  \item{ncmax2}{
integer, set upper bound for \code{k}
}
  \item{ncd}{
integer, set number of conditional densities computed
}
  \item{ngrid}{
integer, set number of grid points for density evaluation
}
  \item{kkzz}{
integer, set number 
}
  \item{idebug}{
integer, number of sweep from which to print debugginh information
}
  \item{qdebug}{
integer 1 or 0 according to whether debugging information is to be printed
}
}
\value{
An object of class \code{nmix}. List with numerous components, including
  \item{post}{posterior distribution of number of components \code{k}}
  \item{pe}{list whose \code{k}'th component is a \code{k} by 3 matrix of estimated posterior means of weights, means and sd's for a mixture with \code{k} components}
  \item{den}{matrix of density estimates for \code{k=1,2,...,6} and overall, preceded by row of abcissae at which they are evaluated - only when \code{out} includes "d"}
  \item{avn}{order-\code{ncmax2} square matrix with \code{(i,j)} entry the posterior expected number of observations allocated to component \code{i} when there are \code{j} components in the mixture - only when \code{out} includes "d"}
  \item{traces}{list of named vectors, traces of selected statistics \code{k}, \code{entropy} (as defined in Green and Richardson, 2001), etc, sub-sampled to every \code{nspace} sweeps}
%% ...
}
\references{
Richardson, S. and Green, P. J. On Bayesian analysis of mixtures with an unknown number of components (with
discussion), J. R. Statist. Soc. B, 1997, 59, 731-792; see also the
correction in J. R. Statist. Soc. B, 1998, 60, 661.

Green, P. J. and Richardson, S. Modelling heterogeneity with and without the Dirichlet process,
Scandinavian Journal of Statistics, 2001, 28, 355-375.
}

\author{
Peter J. Green
}
\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000,out="d")
print(z)
summary(z)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
