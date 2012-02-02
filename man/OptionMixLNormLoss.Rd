\name{OptionMixLNormLoss}
\alias{OptionMixLNormLoss}
\title{
Loss function for fitting n log-normal distributions
}
\description{
This function computes the loss function for fitting n log-normal distributions from option data.
}
\usage{
OptionMixLNormLoss(pars, Probi, k, r, X, F, C, W)
}
\arguments{
  \item{pars}{
nx2 vector,parameters in log-normal distributions in the following order:
(mean_1, mean_2, ..., mean_n, Std_1, Std_2, ..., Std_n)
}
  \item{Probi}{
nx1 vector, weights (probabilites) the different normal distributions
}
  \item{k}{
scalar, time to expiry
}
  \item{r}{
scalar, interest rate per period
}
  \item{X}{
Mx1 vector, option strike prices
}
  \item{F}{
scalar or Mx1 vector, forward prices. If the dimension is Mx1, it contains one forward price for each call price.
}
  \item{C}{
Mx1 vector, call option prices
}
  \item{W}{
(1+M)x(1+M) or (M+M)x(M+M), weighting matrix
}
}
\details{
See Soderlind and Svensson (1997).
}
\value{
\item{loss }{The value of [(model_F-data_F); (model_C-data_C)]'W*[(model_F-data_F);(model_C-data_C)]}, where 'model_' and 'data_' denote, respectively, the fitted and the sample prices. 
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



\examples{
## A simple example:
  
  #strike prices
  X <- matrix( c( 92.00, 94.00, 94.50, 95.00, 95.50, 96.00, 96.50, 97.00,
	97.50, 98.00, 98.50, 99.00, 99.50, 100.0, 100.5, 101.0,
	101.5, 102.0, 102.5, 103.0, 103.5 ) )
	
  #call prices; Mx1 vector
  C <- matrix( c( 5.13, 3.25, 2.83, 2.40, 2.00, 1.64, 1.31, 1.02,
	0.770,  0.570,  0.400,  0.280,  0.190, 0.130, 0.0800, 0.0500,
	0.0400, 0.0300, 0.0200, 0.0100, 0.0100 ) )

  #Forward price, scalar or Mx1 vector
  F <- 97.05                

  #time to expiry in years, scalar
  k <- 48/365               

  #Interest rate: LIFFE margining<->no discounting
  r <- 0                    

  #weight matrix for price errors
  W <- diag( length(F)+length(X) )  

  par0  <- matrix(c( 4.575, 0.03 ))  
  Probi <- 1
  
  results <- optim(par0, OptionMixLNormLoss, gr=NULL, 
		method = c("L-BFGS-B"), lower=0, Probi, k, r, X, F, C, W)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ OptionMixLNormLoss }
