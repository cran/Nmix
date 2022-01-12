\name{galx}
\alias{galx}
\docType{data}
\title{
Galaxy data set
}
\description{
Velocities of 82 distant galaxies, diverging from our own galaxy
}
\usage{data("galx")}
\format{
  The format is:
 num [1:82] 9.17 9.35 9.48 9.56 9.78 ...
}
\source{
Roeder, K. (1990) Density estimation with confidence sets exemplified by superclusters and voids in the
galaxies. J. Am. Statist. Ass., 85, 617-624.
}
\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792.
}
\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000,out="d")
}
\keyword{datasets}
