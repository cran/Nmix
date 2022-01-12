\name{enz}
\alias{enz}
\docType{data}
\title{
Enzyme data set
}
\description{
Enzymatic activity in the blood, for an enzyme involved in the
metabolism of carcinogenic substances, among a group of 245 unrelated individuals.
}
\usage{data("enz")}
\format{
  The format is:
 num [1:245] 0.13 0.08 1.261 0.224 0.132 ...
}
\source{
Bechtel, Y. C., Bonaiti-Pellik, C., Poisson, N., Magnette, J. and Bechtel, P. R. (1993) A population and
family study of N-acetyltransferase using caffeine urinary metabolites. Clin. Pharm. Therp., 54, 134-
141.
}
\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792.
}
\examples{
data(enz)
z<-Nmix('enz',nsweep=5000,nburnin=500,out="d")
}
\keyword{datasets}
