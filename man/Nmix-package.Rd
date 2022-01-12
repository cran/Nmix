\name{Nmix-package}
\alias{Nmix-package}
\docType{package}
\title{
\packageTitle{Nmix}
}
\description{
A program for Bayesian analysis of univariate normal mixtures, implementing the approach of Richardson and Green (1997)
}
\details{
A program for Bayesian analysis of univariate normal mixtures with an
unknown number of components, implementing the approach of Richardson and Green, 
Journal of the Royal Statistical Society, B, 59, 731-792 (1997);
see also the correction in J. R. Statist. Soc. B, 1998, 60, 661).
Computation is by reversible jump Markov chain Monte Carlo; this package is essentially an R interface to the Fortran program
originally written in 1996 for the MCMC sampling, together with some facilities for displaying and summarising the resulting posterior distribution, and reporting the sampler performance.}
\author{
\packageAuthor{Nmix}

Maintainer: \packageMaintainer{Nmix}
}

\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792; see also the
correction in J. R. Statist. Soc. B, 1998, 60, 661)}

\keyword{ package }
\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000,out="d")
print(z)
summary(z)
}