\name{lnacid}
\alias{lnacid}
\docType{data}
\title{
Lake acidity data set
}
\description{
Acidity index measured in a sample of 155 lakes in north-central Wisconsin, on log scale.
}
\usage{data("lnacid")}
\format{
  The format is:
 num [1:155] 2.93 3.91 3.73 3.69 3.82 ...
}
\source{
Crawford, S. L., DeGroot, M. H., Kadane, J. B. and Small, M. J. (1992) Modeling lake chemistry
distributions: approximate Bayesian methods for estimating a finite mixture model. Technometrics,
34, 441-453. 
}
\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792.
}
\examples{
data(lnacid)
z<-Nmix('lnacid',nsweep=10000,nburnin=1000,out="d")
}
\keyword{datasets}
