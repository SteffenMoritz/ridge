\name{linearRidgeGenotypesPredict}
\alias{linearRidgeGenotypesPredict}

\title{
  Predict phenotypes from genome-wide SNP data based on a file of coefficients
}
\description{
  Predict phenotypes from genome-wide SNP data based on a file of
  coefficients. Genotypes and fitted coefficients are provided as
  filenames, allowing the computation of fitted probabilities when SNP
  data are too large to be read into R.
}
\usage{
linearRidgeGenotypesPredict(genotypesfilename, betafilename, phenotypesfilename = NULL,
verbose = FALSE)
}

\arguments{
  \item{genotypesfilename}{
    character string: path to file containing SNP genotypes coded 0, 1,
  2. See \code{Input file formats}.
}
  \item{betafilename}{
    character string: path to file containing fitted coefficients. See \code{Input file formats}.
}
  \item{phenotypesfilename}{
    (optional) character string: path to file in which to write out the
  predicted phenotypes. See \code{Output file formats}. Whether or not this argument
  is supplied, the fitted coefficients are also returned by the function.
}
  \item{verbose}{
    Logical: If \code{TRUE}, additional information is printed to the R
  outupt as the code runs. Defaults to \code{FALSE}.
}
}

\section{Input file formats}{
\describe{
\item{genotypesfilename:}{A header row, plus one row for each
  individual, one SNP per column. The header row contains SNP
  names. SNPs are coded as 0, 1, 2 for minor allele count. Missing
  values are not accommodated. }
\item{betafilename:}{Two columns: First column is SNP names in same order as in \code{genotypesfilename}, second column is fitted coefficients. If the coefficients include an intercept then the first row of \code{betafilename} should contain it with the name Intercept in the first column. An Intercept thus labelled will be used appropriately in predicting the phenotypes. SNP names must match those in \code{genotypesfilename}. 
  The format of \code{betafilename} is
  that of the output of \code{\link{linearRidgeGenotypes}}, meaning
  \code{linearRidgeGenotypesPredict} can be used to predict using
  coefficients fitted using \code{\link{linearRidgeGenotypes}} (see the example).
}
}
}

\section{Output file format}{
Whether or not \code{phenotypesfilename} is provided, predicted phenotypes are returned to the R workshpace. If \code{phenotypesfilename} is provided, predicted phenotypes are written to the file specified (in addition).
\describe{
\item{phenotypesfilename:}{One column, containing predicted phenotypes, one individual per row.}
}
}

\value{
  A vector of fitted values, the same length as the number of
  individuals whose data are in \code{genotypesfilename}. If
  \code{phenotypesfilename} is supplied, the fitted values are also
  written there.
}
\references{
A semi-automatic method to guide the choice of ridge parameter in ridge regression. Cule, E. and De Iorio, M. (2012) arXiv:1205.0686v1 [stat.AP] 
}
\author{
  Erika Cule
}
\seealso{
\code{\link{linearRidgeGenotypes}} for model
  fitting. \code{\link{logisticRidgeGenotypes}} and
  \code{\link{logisticRidgeGenotypesPredict}} for corresponding functions
  to fit and predict on SNP data with binary outcomes.
}
\examples{
\dontrun{
genotypesfile <- system.file("extdata","GenCont_genotypes.txt",package = "ridge")
phenotypesfile <- system.file("extdata","GenCont_phenotypes.txt",package = "ridge")
betafile <- tempfile(pattern = "beta", fileext = ".dat")
beta_linearRidgeGenotypes <- linearRidgeGenotypes(genotypesfilename = genotypesfile,
                                                      phenotypesfilename = phenotypesfile,
                                                      betafilename = betafile)
pred_phen_geno <- linearRidgeGenotypesPredict(genotypesfilename = genotypesfile,
                                                    betafilename = betafile)
## compare to output of linearRidge
data(GenCont) ## Same data as in GenCont_genotypes.txt and GenCont_phenotypes.txt
beta_linearRidge <- linearRidge(Phenotypes ~ ., data = as.data.frame(GenCont))
pred_phen <- predict(beta_linearRidge)
print(cbind(pred_phen_geno, pred_phen))
## Delete the temporary betafile
unlink(betafile)
}
}
