\name{plot.nmix}
\alias{plot.nmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting for Bayesian Inference on Univariate Normal Mixtures
}
\description{
Plotting of various information from \code{nmix} object on current graphics device
}
\usage{
\method{plot}{nmix}(x, which=1, offset=1, nsamp=50, equi=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{nmix} object, as output by \code{Nmix} function
}
  \item{which}{
integer vector, specifying which of several available plots are required, see 'Details' below
}
  \item{offset}{
numeric, vertical displacement between plotted traces, in the case \code{which} contains 2
}
  \item{nsamp}{
integer, number of posterior samples plotted, in the case \code{which} contains 3
}
  \item{equi}{
logical, should thinning of posterior density samples be equi-spaced, rather than randomm, in the case \code{which} contains 3
}
  \item{...}{
additional arguments to \code{\link{plot}}
}
}
\value{
\code{NULL, invisibly}; plot method for class \code{nmix}. Function called for its side effect of plotting selected information from \code{x} on current graphics device
}
\details{
If \code{which} includes 1, a 2-panel plot of which: the first is a histogram of the data, and if \code{z} has a component \code{den}, superimposed plots of the posterior density estimates, conditional on \code{k=1,2,...,6} and unconditionally; and the second a barplot of the estimated posterior distribution of \code{k}.

If \code{which} includes 2, a multiple trace plot of various statistics for a thinned subsample of the MCMC run, after burnin-in. The statistics are the (named) components of the list \code{z$traces}, currently (i) the number of components \code{k}, (ii) the entropy, and (iii) if the deviance trace has been saved (\code{out} includes "d"), the deviance, of the current sample. The traces may be of different lengths, the horizontal scales in the plot are adjusted to span the length of the (post burn-in) MCMC run, regardless of these differences.

If \code{which} includes 3, a thinned sample of size \code{nsamp} from the posterior distribution of the density function, computed from a thinned sample of (weight, mean, sd) values generated in the posterior simulation. 
}
\references{
Richardson and Green (J. R. Statist. Soc. B, 1997, 59, 731-792; see also the
correction in J. R. Statist. Soc. B, 1998, 60, 661)}

\author{
Peter J. Green
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(galx)
z<-Nmix('galx',nsweep=10000,nburnin=1000,out="d")
plot(z,1:2)
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
