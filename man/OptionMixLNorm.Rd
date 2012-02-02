\name{OptionMixLNorm}
\alias{OptionMixLNorm}
\title{
Computes call prices for a vector of option strike prices
}
\description{
This function uses a generic option price formula based on the stochastic discout factor approach: -kr=Em+Var(m)/2. The model accounts for n different states where prices are log-normally distributed. Hence, the pricing formula is based on a mixture of n lognormal distributions. The purpose is compute call prices for M different strike prices. The user should notice that this routine does not use the asset price for the underlying asset. 
}
\usage{
OptionMixLNorm(k, r, sbar_ms, SigmaSS, X, Probi)
}
\arguments{
  \item{k}{
scalar, time to expiry
}
  \item{r}{
scalar, interest rate per period
}
  \item{sbar_ms}{
nx1 vector, E(t)lnS(t+k) + Cov(t)[lnM(t+k),lnS(t+k)]
}
  \item{SigmaSS}{
nx1 vector, Var(t)[lnS(t+k)]
}
  \item{X}{
Mx1 vector, option strike prices
}
  \item{Probi}{
nx1 vector of probabilities of different states
}
}
\details{
See Soderlind and Svensson (1997).
}
\value{
\item{C }{Mx1 vector of call prices at M strike prices}
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
  
  #strike prices
  X <- matrix( c( 92.00, 94.00, 94.50, 95.00, 95.50, 96.00, 96.50, 97.00,
	97.50, 98.00, 98.50, 99.00, 99.50, 100.0, 100.5, 101.0,
	101.5, 102.0, 102.5, 103.0, 103.5 ) )

  #time to expiry in years, scalar
  k <- 48/365               

  #Interest rate: LIFFE margining<->no discounting
  r <- 0  

  Probi 		<- 1
  pars  		<- matrix(c( 4.575, 0.03 ))
  n       		<- floor( (length(pars))/2 )      
  sbar_ms <- matrix(c( pars[1:n] )) 

  #estimating standard deviation
  g1			<- n+1
  g2			<- 2*n
  SigmaSS <- matrix(c( pars[g1:g2]^2 ))
  
  y <- OptionMixLNorm(k, r, sbar_ms, SigmaSS, X, Probi)  
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ OptionMixLNorm }

