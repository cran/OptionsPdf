

erf <- function(x){
return(2*pnorm(x*sqrt(2))-1)
}




OptionMixLNormLoss <- function(pars, Probi, k, r, X, F, C, W) {


#  OptionMixLNormLoss    Loss function for fitting n log-normal distributions from 
#                        option data.
#
#  Input:     pars        nx2 vector of parameters in log-normal
#                         distributions in the following order:
#                         mean_1,mean_2,...,mean_n,
#                         Std_1 ,Std_2 ,...,Std_n,
#             Probi       nx1 vector of weights (probabilites) the different
#                         normal distributions
#             k           scalar, time to expiry
#             r           scalar, interest rate per period
#             X           Mx1 vector, strike prices
#             F           scalar or Mx1 vector, forward price. If Mx1,
#                         then it is one forward price for each call price.
#             C           Mx1 vector, call option prices
#             W           (1+M)x(1+M) or (M+M)x(M+M) weighting matrix
#
#  Output:    loss	  the value of [(F-g_F);(C-g_C)]'g_W*[(F-g_F);(C-g_C)]


#states
n <- floor( (length(pars))/2 )
               
sbar_ms <- matrix(c( pars[1:n] )) 


#estimating standard deviation
g1=n+1
g2=2*n
SigmaSS <- matrix(c( pars[g1:g2]^2 ))


#implied call price
Chat <- OptionMixLNorm(k, r, sbar_ms, SigmaSS, X, Probi)   


#implied forward rate
Fhat <- sum(c( Probi * exp(sbar_ms + SigmaSS/2) ) )        

dev  <- matrix( c(Fhat-F, Chat-C) , byrow=TRUE)

loss <- t(dev) %*% W %*% dev * 10000

return(loss)

}






OptionMixLNorm <- function(k, r, sbar_ms, SigmaSS, X, Probi) {


# OptionMixLNorm     Generic option price formula, using -kr=Em+Var(m)/2 and n
#                    different states (mixture of n lognormal distributions).
#                    Does not use asset pricing equation for underlying asset.
#                    Calculates call prices for M different strike prices.
#
#  Inputs:    k             scalar, periods to expiry date  (could be <1 );
#             r             scalar, interest rate per period
#             sbar_ms       nx1 vector of E(t)lnS(t+k) +
#                           Cov(t)[lnM(t+k),lnS(t+k)]
#             SigmaSS       nx1 vector of Var(t)[lnS(t+k)]
#             X             Mx1 vector of strike prices
#             Probi         nx1 vector of probabilities of different states
#
#  Output:    C             Mx1 vector of call prices at M strike prices
#
#  IMPORTANT: requires the library Matlab


library(matlab)

n <- length(Probi) 
M <- matrix(c( length(X) ))

sbar_ms <- repmat( sbar_ms, 1, M )     
  
SigmaSS <- repmat( SigmaSS, 1, M )
Probi   <- repmat( Probi,   1, M )

X       <- repmat( c(t(X)), n, 1)

d1 <- ( sbar_ms + SigmaSS - log(X) )/sqrt(SigmaSS)

d2 <- ( sbar_ms -           log(X) )/sqrt(SigmaSS)

#std normal cdf of d1 and d2
Ncdfd1 <- 0.5 + 0.5*erf( d1/sqrt(2) )         
Ncdfd2 <- 0.5 + 0.5*erf( d2/sqrt(2) )


C1 <- exp( -k*r )*exp(sbar_ms + SigmaSS/2) * Ncdfd1 - exp( -k*r )*X * Ncdfd2 

#expectation over states -> Mx1 vector
Chat <- matrix(c(colSums(Probi*C1)))  


}





MixedNormalPdf <- function(x, mu, s2, weight) {


# MixedNormalPdf: Calculating univariate pdf of mixture of normal distributions.
#
#  Input:     x        TxK matrix
#             mu       nx1 or 1xn vector, means of n different N(mu,s2) distributions
#             s2       nx1 or 1xn vector, variances of n different N(mu,s2) distributions
#             weight   nx1 or 1xn vector, weight of n different N(mu,s2) distributions.
#                      Should sum to unity
#
#  Output:    pdf      TxK matrix, pdf of each of the columns of x
#             pdf_i    TxKxn matrix or Txn (if K=1), pdf_i(:,:,i) is the pdf of mixture component i
#
#
# IMPORTANT: requires the library Matlab


library(matlab)


#weights should sum to unity
if (abs(sum(weight)-1) > 1e-9) {     
	warning("Weights do not sum to unity")
end
}


T <- nrow(x) 
K <- ncol(x) 


#no. of mixed normal distributions
n <- length(weight)          


pdf <- c(0);

i <- 1

#loop over components (normal distributions)
for (i in 1:n) {     
	#standardized x       
	xStd <- matrix(c( (x-mu[i])/sqrt(s2[i]) ))         
	pdf  <- matrix(c( pdf + weight[i] * exp( -xStd^2/2 )/sqrt(2*pi*s2[i]) ))
}


pdf_i <- repmat(NaN, T, K, n);

i <- 1

#loop over components (normal distributions)
for (i in 1:n) {  
        #standardized x     
	xStd         <- matrix(c( (x-mu[i])/sqrt(s2[i]) ))         
	pdf_i[ , ,i] <- t(c( exp( -(xStd^2)/2 )/sqrt(2*pi*s2[i]) ))
}


#TxKxn -> Txn if K=1
pdf_i <- drop(pdf_i) 


return(list(pdf=pdf,pdf_i=pdf_i))

}


