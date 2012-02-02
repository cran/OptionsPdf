\name{OptionsPdf-package}
\alias{OptionsPdf-package}
\alias{OptionsPdf}
\docType{package}
\title{
This package estimates the implied-probability density function from interest-rate option price data.
}
\description{
This package provides the estimation of implied probability density functions (pdf) from data on interest-rate option prices. The pdf's take the form of a mix of lognormal distributions. The methods used in the routines are introduced by Soderlind and Svensson (1997) 
The user should kindly keep in mind that this package is the R port of the Matlab functions originally distributed by P. Soderlind.
}
\details{
\tabular{ll}{
Package: \tab OptionsPdf\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-02\cr
License: \tab GPL-3\cr
}
The key function is OptionMixLNormLoss.R,. This generates the loss function value from fitting an n-dimensional log-normal distribution to option price data. The optimizing parameters consist of means and standard deviations for a mixture of normal distributions. These parameters minimize the weighted squared difference between model-consistent and observed prices of options. The loss function takes into account both option call prices and the associated forward prices. The function OptionMixLNorm uses a formula that acconts for n different states of the world to compute the call prices from a vector of strike prices. The function MixedNormalPdf delivers the univariate probability density functions from a mixture of normal distributions and, thus, allows the user to plot the density functions implied by the option market data. 
}

\author{
Paolo Zagaglia

Maintainer: Paolo Zagaglia <paolo.zagaglia@gmail.com>
}

\references{
Soderlind, P., and L. E. O. Svensson (1997), "New Techniques to Extract Market Expectations from Financial Instruments", Journal of 	
Monetary Economics, 40, 383-429.
}

\examples{
#strike prices
X <- matrix( c(92.00, 94.00, 94.50, 95.00, 95.50, 96.00, 96.50, 97.00,
	97.50, 98.00, 98.50, 99.00, 99.50, 100.0, 100.5, 101.0,
	101.5, 102.0, 102.5, 103.0, 103.5 ) )
	
#call prices; Mx1 vector
C <- matrix( c( 5.13,   3.25,   2.83,   2.40,  2.00,  1.64,   1.31,   1.02,
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


#  Two states: solving for 2 means, 2 standard deviations, and
#  one probability of the first normal distribution.
#
#  Estimation: (a) Assuming a probability and then optimizing wrt.
#                     mean1, mean2, std1, and std2.
#                     (b) (a) is repeated for (Nprob) different values of
#                     the probability between 0.01 and 0.49. 
#
#  Uses forward rate and option data.
  
#initial parameter guess
par0 <- c(4.58, 4.58, 0.025, 0.005)  

#Number of different probabilites in [0.01,0.49] to try
#NProb x 5 matrix to store loop results in
NProb <- 43          
                     
#probabilities to try             
Prob   <- matrix(c( linspace(0.03,0.47,NProb) ))   
Par2M  <- repmat( NaN,NProb,4)

#NProb x 1 vector to store loss fn values
Loss2M <- repmat( NaN,NProb,1 )   

i <- 1;                              
  #looping over different probabilities
  for (i in 1:NProb) {

	#[Prob(i);1-Prob(i)]
	Probi <- matrix(c( Prob[i], 1-Prob[i] ) )
  
	results2 <- optim(par0, OptionMixLNormLoss, 
			gr=NULL, Probi, k, r, X, F, C, W)
 
	Par2M[i,] <- t(c( results2$par )) 
	Loss2M[i] <- t(c( results2$value ))
  }

  #index of iteration with minimum loss
  MinLossIndex <- find( Loss2M==min(Loss2M) )

  #Parameters at minimum loss 
  par2 <- matrix(c(Par2M[MinLossIndex,], Prob[MinLossIndex]))

  S <- linspace(4.47, 4.68, 100)

  pdf2 <- MixedNormalPdf( matrix(S), matrix(c(par2[1:2])), 
		matrix(c(par2[3:4]^2)), matrix(c(par2[5],1-par2[5])) )

  plot( S, pdf2$pdf, type="b", col="red",  xlab="log Bund price", 
		main="Pdf from LIFFE Bunds option data, 6 April 1994 - 2")   

}
