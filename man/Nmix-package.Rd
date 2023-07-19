\name{Nmix-package}
\alias{Nmix-package}
\docType{package}
\title{
\packageTitle{Nmix}
}
\description{
A program for Bayesian analysis of univariate normal mixtures, implementing the approach of Richardson and Green (1997) <doi:10.1111/1467-9868.00095>
}
\details{
A program for Bayesian analysis of univariate normal mixtures with an
unknown number of components, implementing the approach of Richardson and Green, 
Journal of the Royal Statistical Society, B, 59, 731-792 (1997);
see also the correction in J. R. Statist. Soc. B, 1998, 60, 661).
Computation is by reversible jump Markov chain Monte Carlo; this package is essentially an R interface to the Fortran program
originally written in 1996 for the MCMC sampling, together with some facilities for displaying and summarising the resulting posterior distribution, and reporting the sampler performance.}
\author{
Maintainer: \packageMaintainer{Nmix}
}

\references{
Richardson and Green (1997) <doi:10.1111/1467-9868.00095> (J. R. Statist. Soc. B, 59, 731-792; see also the
correction in <doi:10.1111/1467-9868.00146>, J. R. Statist. Soc. B, 1998, 60, 661).

The author is grateful to Peter Soerensen for providing the interface to the C i/o routines used here.
}
\keyword{ package }
\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000)
print(z)
summary(z)
}
