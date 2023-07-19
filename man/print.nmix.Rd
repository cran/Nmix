\name{print.nmix}
\alias{print.nmix}
\title{Printing for Bayesian Inference on Univariate Normal Mixtures}
\description{Printing of various information about \code{nmix} object on current output}
\usage{\method{print}{nmix}(x, ...)}
\arguments{
  \item{x}{\code{nmix} object, as output by \code{Nmix} function}
  \item{...}{additional arguments to \code{\link{print}}}
}
\value{
\code{x}, invisibly; print method for class \code{nmix}. Function called for its side effect of printing selected information from \code{x}
}
\details{
Currently the information printed consists of the estimated posterior for \code{k} and basic parameters 
of the MCMC simulation: number of sweeps, length of burnin, random number seed to replicate the run, and confirmation of which MCMC moves used (latters s,w,p,a,h,b standing for split/merge, weights, parameters, allocations, hyperparameters and birth/death).
}
\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792; see also the
correction in J. R. Statist. Soc. B, 1998, 60, 661)
}
\author{
Peter J. Green
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000,out="d")
z
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
