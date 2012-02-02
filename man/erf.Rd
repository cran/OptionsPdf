\name{erf}
\alias{erf}
\title{
Error function
}
\description{
This routine returns the value of the error function. This function is twice the integral of the Gaussian distribution with 0 mean and variance of 1/2.
}
\usage{
erf(x)
}
\arguments{
  \item{x}{
Real-valued vector
}
}
\value{
\item{y }{Real-valued vector}
}


\author{
Paolo Zagaglia, paolo.zagaglia@gmail.com
}


\examples{
## A simple example:
x<-c(0.4)
y<-erf(x)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ erf }
