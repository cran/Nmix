\name{sdrni}
\alias{sdrni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random number initialiser, allowing retrospective replication
}
\description{
Front-end to standard R random number seed setter, allowing retrospective replication
}
\usage{
sdrni(seed) 
}
\arguments{
  \item{seed}{
non-negative integer random number seed, often 0 for absolute re-initialisation as with \code{set.seed(NULL)}
}
}
\value{
\code{seed} if input value is positive, otherwise the value that if used in a subsequent call will deliver exactly the same random numbers
}
\details{
Using \code{sdrni} to initialise random number stream allows a decision to repeat a simulation exactly, presumably with additional outputs, need only be made after seeing results; see Examples
}
\author{
Peter J. Green
}

\examples{
sdrni(0)
runif(5)
keep<-sdrni(0)
runif(5)
sdrni(keep)
runif(5)
}

