\name{kripp.boot}
\alias{kripp.boot}

\title{
Bootstrap Krippendorff's alpha for intercoder reliability data
}
\usage{
kripp.boot(x, raters='rows', probs=c(.025,.975), iter=100, method='nominal')
}
\arguments{
\item{x}{Dataframe or matrix containing classifications or scores}
\item{raters}{Is each rater's data on the column or row?}
\item{probs}{2-item vector of probabilities with values in [0,1]}
\item{iter}{The number of bootstrap resamplings to be performed}
\item{method}{Data level of x}

\details{
Convenience wrapper for the kripp.alpha function in the irr package. Performs bootstrapping (resampling with replacement) on intercoder reliability measures, returns confidence intervals of Krippendorff's alpha. 
}

\references{
Krippendorff, K. (2003). Content analysis: An introduction to its methodology. Beverly Hills, CA: Sage.
}

\author{
Mike Gruszczynski, <mikewgruz@gmail.com>
}
\examples{
# coder data from Krippendorff
nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)

# each coder is on a row, so assume default rater value
bootstrapped.alpha <- kripp.boot(nmm, iter=1000)
}


