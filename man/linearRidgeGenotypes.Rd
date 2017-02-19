\name{linearRidgeGenotypes}
\alias{linearRidgeGenotypes}

\title{
  Fits linear ridge regression models for genome-wide SNP data.
}
\description{
Fits linear ridge regression models for genome-wide SNP data. The SNP
  genotypes are not read into R but file names are passed the code
  directly, enabling the analysis of genome-wide scale SNP data sets.
}
\usage{
linearRidgeGenotypes(genotypesfilename, phenotypesfilename, lambda = -1, 
thinfilename = NULL, betafilename = NULL, approxfilename = NULL,
permfilename = NULL, intercept =
TRUE, verbose = FALSE)
}

\arguments{
  \item{genotypesfilename}{
    character string: path to file containing SNP genotypes coded 0, 1,
  2. See \code{Input file formats}.
  }
  \item{phenotypesfilename}{
    character string: path to file containing phenotypes. See \code{Input file formats}.
  }
  \item{lambda}{
    (optional) shrinkage parameter. If not provided, the default denotes
  automatic choice of the shrinkage parameter using the method of Cule &
  De Iorio (2012). 
}
  \item{thinfilename}{
(optional) character string: path to file containing three columns: SNP name, chromosme and SNP psotion. See \code{Input file formats}. (See \code{details}.)
}
  \item{betafilename}{
(optional) character string: path to file where the output will be written. See \code{Output file formats}.
}
\item{approxfilename}{
(optional) character string: path to fine where the approximate test p-values will be written. 
Approximate p-values are not computed unless this argument is given. Approximate p-values
  are computed using the method of Cule et al (2011). See \code{Output file formats}.
}
\item{permfilename}{
  (optional) character string: path to file where the permutation test
  p-values will be written.
  Permutation test p-values are not computed unless this argument is
  given. (See warning). See \code{Output file formats}.
}
\item{intercept}{
  Logical: Should the ridge regression model be fitted with an
  intercept? (Defaults to \code{TRUE})
}
\item{verbose}{
  Logical: If \code{TRUE}, additional information is printed to the R
  output as the code runs. Defaults to \code{FALSE}.
}
}

\section{Input file formats}{
\describe{
\item{genotypesfilename:}{A header row, plus one row for each
  individual, one SNP per column. The header row contains SNP
  names. SNPs are coded as 0, 1, 2 for minor allele count. Missing
  values are not accommodated. Invariant SNPs in the data cause an
  error, please remove these from the file before calling the function.}
\item{phenofilename:}{A single column of phenotypes with the individuals in the same order as those in the file \code{genotypesfilename}.}
\item{thin:}{(optional) Three columns and the same number of rows as there are SNPs in the file \code{genotypesfilename}, one row per SNP. First column: SNP names (must match names in \code{genotypesfilename}); second column: chromosome;  third column: SNP position in BP.}
}
}

\section{Output file formats}{
All output files are optional. Whether or not \code{betafilename} is provided, fitted coefficients are returned to the R workshpace. If \code{betafilename} is provided, fitted coefficients are written to the file specified (in addition).
\describe{
\item{betafilename:}{Two columns: First column is SNP names in same order as in \code{genotypesfilename}, second column is fitted coefficients. If \code{intercept = TRUE} (the default) then the first row is the fitted intercept (with the name Intercept in the first column).}
\item{approxfilename:}{Two columns: First column is SNP names in same order as in \code{genotypesfilename}, second column is approximate p-values.}
\item{permfilename:}{Two columns: First column is SNP names in same order as in \code{genotypesfilename}, second column is permutation p-values.}
}
}

\details{
  If a file \code{thin} is supplied, and the shrinkage parameter
  \code{lambda} is being computed automatically based on the data, then
  this file is used to thin the SNP data by SNP position. If this file
  is not supplied, SNPs are thinned automatically based on number of SNPs.
}

\value{
  The vector of fitted ridge regression coefficients.
  If \code{betafilename} is given, the fitted coefficients are written to this
  file as well as being returned.
  If \code{approxfilename} and/or \code{permfilename} are given, results of approximate
  test p-values and/or permutation test p-values are written to the files
  given in their arguments.
}

\section{Warning }{
  When data are large, the permutation test p-values
  may take a very long time to compute. It is recommended not to request
  permutation test p-values (using the argument \code{permfilename})
  when data are large.}

\references{
  Significance testing in ridge regression for genetic data. Cule, E. et
  al (2011) BMC Bioinformatics, 12:372
  A semi-automatic method to guide the choice of ridge parameter in
  ridge regression. Cule, E. and De Iorio, M. (2012) arXiv:1205.0686v1 [stat.AP]
}

\author{
  Erika Cule
}

\seealso{
  \code{\link{linearRidge}} for fitting linear ridge regression models
  when the data are small enough to be read into R.
  \code{\link{logisticRidge}} and \code{\link{logisticRidgeGenotypes}} for fitting logistic ridge
  regression models.
}

\examples{
\dontrun{
    genotypesfile <- system.file("extdata","GenCont_genotypes.txt",package = "ridge")
    phenotypesfile <- system.file("extdata","GenCont_phenotypes.txt",package = "ridge")
    beta_linearRidgeGenotypes <- linearRidgeGenotypes(genotypesfilename = genotypesfile,
phenotypesfilename = phenotypesfile)
    ## compare to output of linearRidge
    data(GenCont) ## Same data as in GenCont_genotypes.txt and GenCont_phenotypes.txt
    beta_linearRidge <- linearRidge(Phenotypes ~ ., data = as.data.frame(GenCont))
    cbind(round(coef(beta_linearRidge), 6), beta_linearRidgeGenotypes)
}
  }
