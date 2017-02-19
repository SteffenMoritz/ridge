\name{pvals}
\alias{pvals}
\alias{pvals.ridgeLinear}
\alias{pvals.ridgeLogistic}
\alias{plot.pvalsRidgeLinear}
\alias{plot.pvalsRidgeLogistic}
\alias{print.pvalsRidgeLinear}
\alias{print.pvalsRidgeLogistic}
\title{
Compute p-values for ridgeLinear and ridgeLogistic models
}
\description{
Functions for computing, printing and plotting p-values for ridgeLinear
and ridgeLogistic models. The p-values are computed using the
significance test of Cule et al (2011). 
}
\usage{
pvals(x, ...)

\method{pvals}{ridgeLinear}(x, ...)

\method{pvals}{ridgeLogistic}(x, ...)

\method{print}{pvalsRidgeLinear}(x, digits = max(3, getOption("digits") - 3),
signif.stars = getOption("show.signif.stars"), all.coef = FALSE, ...) 

\method{print}{pvalsRidgeLogistic}(x, digits = max(3, getOption("digits") - 3),
signif.stars = getOption("show.signif.stars"), all.coef = FALSE, ...)

\method{plot}{pvalsRidgeLinear}(x, y = NULL, ...)

\method{plot}{pvalsRidgeLogistic}(x, y = NULL, ...)

}
\arguments{

  \item{x}{For the pvals methods, an object of class "ridgeLinear" or "ridgeLogistic", typically from a call to "linearRidge" or 
    "logisticRidge". For the print and plot methods, an object of class
    "pvalsRidgeLinear" or "pvalsRidgeLogistic", typically from a
    call to "pvals".}
  \item{digits}{minimum number of significant digits to be used for most numbers}
  \item{signif.stars}{logical; if \code{TRUE}, P-values are additionally encoded
    visually as \code{significance stars} in order to help scanning of
    long coefficient tables.  It defaults to the
    \code{show.signif.stars} slot of \code{options}.
  }
  \item{all.coef}{Logical. Should p-values for all the ridge
    regression parameters be printed, or only the one from the ridge
    parameter chosen using the method of Cule et al (2012)}
  \item{y}{Dummy argument for compatibility with the default \code{plot} method. Ignored.}
  \item{\dots}{further arguments to be passed to or from other methods}

}
\details{
	Standard errors, test statistics and p-values are computed using coefficients and data on the scale that was 
	used to fit them. If the coefficients were standardized before the model was fitted, then the p-values 
	relate to the scaled data.
}
\value{
  For the pvals methods, an object of class "pvalsRidgeLinear" or "pvalsRidgeLogistic" 
  which is a list with elements
  \item{coef}{The (scaled) regression coefficients}
  \item{se}{The standard errors of the regression coefficients}
  \item{tstat}{The test statistic of the regression coefficients}
  \item{pval}{The p-values of the regression coefficients}
  \item{isScaled}{Were the data scaled before the regression
    coefficients were fitted?}
  For the print methods, the argument \code{x} is returned invisibly.
}
\references{
Significance testing in ridge regression for genetic data. Cule, E. et al (2011) BMC Bioinformatics, 12:372
}
\author{
Erika Cule
}
\seealso{
\code{linearRidge}, \code{logisticRidge}
}
\examples{
data(GenBin)
mod <- logisticRidge(Phenotypes ~ ., data = as.data.frame(GenBin))
pvalsMod <- pvals(mod)
print(pvalsMod)
print(pvalsMod, all.coef = TRUE)
plot(pvalsMod)
}
