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
\method{plot}{nmix}(x, which=1:5, separate=FALSE, plugin=FALSE, offset=1, nsamp=50, 
	equi=TRUE, allsort=TRUE, trued=NULL, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{nmix} object, as output by \code{Nmix} function
}
  \item{which}{
integer vector, specifying which of several available plots are required, see 'Details' below
}
  \item{separate}{
logical, if \code{TRUE} opens a fresh default device for each new plot, otherwise prompts before overwriting a previous plot
}
  \item{plugin}{
logical, should plug-in estimator of density, computed from posterior means of parameters, be superimposed on density plot in darkgreen, in the case \code{which} contains 1?
}
  \item{offset}{t
numeric, vertical displacement between plotted traces, in the case \code{which} contains 2.
}
  \item{nsamp}{
integer, number of posterior samples plotted, in the case \code{which} contains 3.
}
  \item{equi}{
logical, should thinning of posterior density samples be equi-spaced, rather than random, in the case \code{which} contains 3
}
  \item{allsort}{
logical, should observations be sorted before making posterior clusters plot, in the case \code{which} contains 4
}
  \item{trued}{
vectorised function defining a probability density function to be superimposed in blue on density plots, in the cases \code{which} contains 1 and/or 3
}
  \item{...}{
additional arguments to \code{\link{plot}}
}
}
\value{
\code{NULL, invisibly}; plot method for class \code{nmix}. Function called for its side effect of plotting selected information about the fitted posterior distribution and sampler performance from \code{x} on the current graphics device
}
\details{
If \code{which} includes 1, a 2-panel plot of which: the first is a histogram of the data, and if \code{z} has a component \code{den} (\code{Nmix} output option D), superimposed plots of the posterior density estimates, conditional on \code{k=1,2,...,6} and unconditionally (in red); and the second a barplot of the estimated posterior distribution of \code{k}.

If \code{which} includes 2, a multiple trace plot of various statistics for a thinned subsample of the MCMC run, after burn-in. The statistics are the (named) components of the list \code{z$traces} that are numerical vectors, some or all of (i) the number of components \code{k} (\code{Nmix} output option k), (ii) the entropy (\code{Nmix} output option e), and (iii) the deviance (\code{Nmix} output option d), of the current sample. The traces may be of different lengths, the horizontal scales in the plot are adjusted to span the length of the (post burn-in) MCMC run, regardless of these differences.

If \code{which} includes 3 (and \code{Nmix} output option p is present), a thinned sample of size \code{nsamp} from the posterior distribution of the density function, computed from a thinned sample of (weight, mean, sd) values generated in the posterior simulation. 

If \code{which} includes 4 (and \code{Nmix} output option a is present), an image plot showing the posterior probabilities that the corresponding observations are in the same mixture component. Not recommended for large \code{n}.

If \code{which} includes 5 (and \code{Nmix} output option C is present), a 4-panel plot displaying Bayersian classifications based on the fitted model, analogous to Fig. 10 in the Richardson and Green paper. The 4 panels corresponding to conditioning on the 4 values of \code{k} most probable according to the posterior (among those for which the necessary posterior sample information has been recorded (see argument \code{k1k2} of \code{Nmix}), and excepting \code{k=1}).
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
