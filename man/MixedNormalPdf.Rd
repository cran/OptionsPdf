\name{MixedNormalPdf}
\alias{MixedNormalPdf}
\title{
Univariate probability density function of mixture of normal distributions
}
\description{
Calculates the univariate pdf of a mixture of normal distributions
}
\usage{
MixedNormalPdf(x, mu, s2, weight)
}
\arguments{
  \item{x}{
TxK matrix, underlying asset prices (e.g., bund prices in the case of bund options)
}
  \item{mu}{
nx1 or 1xn vector, means of n different N(mu,s2) distributions
}
  \item{s2}{
nx1 or 1xn vector, variances of n different N(mu,s2) distributions
}
  \item{weight}{
nx1 or 1xn vector, weight of n different N(mu,s2) distributions. The vector elements should sum to unity.
}
}
\details{
See Soderlind and Svensson (1997).
}
\value{
\item{pdf }{TxK matrix, pdf of each of the columns of x}
\item{pdf_i }{TxKxn matrix or Txn (if K=1), pdf_i[,i] is the pdf of mixture component i}
}
\references{
P. Soderlind and L. E. O. Svensson (1997), "New Techniques to Extract Market Expectations from Financial Instruments", Journal of Monetary Economics
}
\author{
Paolo Zagaglia, paolo.zagaglia@gmail.com
}
\note{
This is the R port of an analogous Matlab function distributed by P. Soderlind.
}


\section{Requested packages }{
This function requires the package matlab. 
}


\examples{
## A simple example:
  
  par0 <- c(4.58, 4.58, 0.025, 0.005)
  
  %vector of prices of underlying asset
  S    <- linspace(4.47, 4.68, 100)
 
  pdf2 <- MixedNormalPdf( matrix(S), matrix(par0[1:2]), 
			matrix(par0[3:4]^2), matrix(c(0.6,0.4)) )
}
  
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MixedNormalPdf }
