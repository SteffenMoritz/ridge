\name{computeRidgeLogistic}
\alias{computeRidgeLogistic}
\alias{updateBeta}
\alias{objectiveFunction}
\title{
Internal functions for logistic ridge regression.
}
\description{
Internal functions for logisitc ridge regression.
}
\usage{
computeRidgeLogistic(X, y, k, intercept = TRUE, doff = FALSE)
updateBeta(B, X, y, k, intercept = TRUE, doff = FALSE) 
objectiveFunction(B, X, y, k, intercept = TRUE)
}
\arguments{
  \item{X}{
	Matrix of predictors.
}
  \item{y}{
	vector of outcomes.
}
  \item{k}{
	ridge regression parameter.
}
  \item{intercept}{
	does the model have an intercept?
}
  \item{doff}{
	should degrees of freedom of the model be computed?
}
}
\details{
	These functions are called in the function \code{logisticRidge}.
	They are not for calling directly by the user.
}
\value{
	\code{computeRidgeLogistic} returns the fitted logistic ridge regression coefficients. If \code{doff = TRUE} 
	it also returns the degrees of freedom of the model and the degrees of freedom for variance.
	\code{updateBeta} returns the fitted coefficients after one iteration of the Newton-Raphson algorithm. If 
	\code{doff = TRUE}, it also returns the penalty matrix and weights matrix used to compute the degrees of 
	freedom. 
	\code{objectiveFunction} returns the objective function for the current iteration of the Newton-Raphson 
	algorithm. 
}
\references{
A semi-automatic method to guide the choice of ridge parameter in ridge regression. Cule, E. and De Iorio, M. (2012) arXiv:1205.0686v1 [stat.AP] 
}
\author{
Erika Cule
}
\note{
These functions are not to be called directly by the user. They should be called via \code{logisticRidge}.
}

\seealso{
\code{\link{logisticRidge}}
}
\keyword{internal}
