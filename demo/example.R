#----------------------------------------------------------------------------------------------------
#
#  Example program for estimating implied distribution from option prices with the 
#  method proposed by Paul Söderlind and Lars E. O. Svensson (1997), 
#  "New Techniques to Extract Market Expectations from Financial Instruments", 
#  Journal of Monetary Economics.
#
#  This example uses LIFFE options on Bunds futures traded on April 6, 1994 
#  as reported in the original Matlab routines of Paul Söderlind
#
#
#  IMPORTANT: uses the Matlab package
#
#----------------------------------------------------------------------------------------------------


library(matlab)


#----------------------------------------------------------------------------------------------------
#
#  LIFFE Bunds option data, trade date April 6, 1994
#

#strike prices; Mx1 vector
X <- matrix( c( 92.00,  94.00,  94.50,  95.00,  95.50,  96.00,  96.50,  97.00,
	97.50,  98.00,  98.50,  99.00,  99.50,  100.0,  100.5,  101.0,
	101.5,  102.0,  102.5,  103.0,  103.5 ) )

#call prices; Mx1 vector
C <- matrix( c( 5.13,    3.25,    2.83,    2.40,    2.00,    1.64,   1.31,    1.02,
	0.770,   0.570,   0.400,   0.280,   0.190,   0.130,  0.0800,  0.0500,
	0.0400,  0.0300,  0.0200,  0.0100,  0.0100 ) )

#Forward price, scalar or Mx1 vector
F <- 97.05                

#time to expiry in years, scalar
k <- 48/365               

#Interest rate: LIFFE margining<->no discounting
r <- 0                    

#weight matrix for price errors
W <- diag( length(F)+length(X) )  



#----------------------------------------------------------------------------------------------------
#  One state: solving for mean and standard deviation of distribution.
#  Uses: forward rate and option data.

#Initial parameter guess: mean, Std
par0  <- matrix(c( 4.575, 0.03 ))         
Probi <- 1 
#----------------------------------------------------------------------------------------------------

results <- optim(par0, OptionMixLNormLoss, gr=NULL, Probi, k, r, X, F, C, W)



#----------------------------------------------------------------------------------------------------
#  Two states: solving for means and standard deviations of distribution.

#Initial parameter guess
par0 <- matrix(c( 4.58, 4.58, 0.025, 0.005 ))

NProb <- 43          
                     
#probabilities to try             
Prob   <- matrix(c( linspace(0.03,0.47,NProb) ))      

i <- 1
Probi <- matrix(c( Prob[i], 1-Prob[i] ) )
#----------------------------------------------------------------------------------------------------

results1 <- optim(par0, OptionMixLNormLoss, gr=NULL, Probi, k, r, X, F, C, W)





#----------------------------------------------------------------------------------------------------
#
#  Two states: solving for 2 means, 2 standard deviations, and
#  one probability of the first normal distribution.
#
#  Estimation:  (a) Assuming a probability and then optimizing wrt.
#                   mean1, mean2, std1, and std2.
#               (b) (a) is repeated for (Nprob) different values of
#                   the probability between 0.01 and 0.49. 
#
#  Uses forward rate and option data.
#
#----------------------------------------------------------------------------------------------------

  
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
  
	results2 <- optim(par0, OptionMixLNormLoss, gr=NULL, Probi, k, r, X, F, C, W)
 
	Par2M[i,] <- t(c( results2$par )) 
	Loss2M[i] <- t(c( results2$value ))
}


#index of iteration with minimum loss
MinLossIndex <- find( Loss2M==min(Loss2M) )

#Parameters at minimum loss
par2 <- matrix(c(Par2M[MinLossIndex,], Prob[MinLossIndex]))

#minimum loss     
loss2 <- Loss2M[MinLossIndex,]  
 

print('Estimated parameters for n=2 (from loop): ')
print(par2)
print('Loss function value (from loop): ')
print(loss2)


#----------------------------------------------------------------------------------------------------
#plotting some results

%vector of Bund option prices
S    <- linspace(4.47, 4.68, 100)

pdf2 <- MixedNormalPdf( matrix(S), matrix(c(par2[1:2])), matrix(c(par2[3:4]^2)), matrix(c(par2[5],1-par2[5])) )

plot( S, pdf2$pdf, type="b", col="red",  xlab="log Bund price", main="Mixture of normals: Pdf from LIFFE Bunds option data") 
#----------------------------------------------------------------------------------------------------




