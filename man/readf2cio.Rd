\name{readf2cio}
\alias{readf2cio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading binary file of structured binary numerical data
}
\description{
Reading binary file of structured binary numerical data, for use in reading into R numerical data written from Fortran
}
\usage{
readf2cio(fn,imax=Inf,verbose=FALSE) 
}
\arguments{
  \item{fn}{
character variable, path to file to be read.
}
  \item{imax}{
maximum number of list components to be read.
}
  \item{verbose}{
boolean, should the reading be reported?
}
}
\value{
numeric list, vector or matrix according to layout of information in \code{fn}; see Details. 
}
\details{
The function is designed to expedite the transfer of possibly large quantities of numeric information, in binary form, written, typically incrementally, in a Fortran routine called from R, without using the arguments to the function.

Assumed format for binary files holding lists, matrices or vectors of numeric data: \cr
writable from Fortran via f2cio interface, readable in R using readBin \cr
file structure supported: binary file, with integer(4), real(4) or double(8) data \cr
first record: list: 0 0 \cr
          matrix or vector: nc mode (mode = 1, 2 or 3 for integer(4), real(4) or double(8)) \cr
succeeding records, one per component of list or row of matrix:  \cr
          list: number of items, mode as integers, followed by data for this component
          (note that modes can differ between but not within components) \cr
          matrix or vector: data for this row \cr
          one-column matrices are delivered as vectors \cr
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
