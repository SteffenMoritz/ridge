\name{GenCont}
\alias{GenCont}
\alias{GenCont_genotypes}
\alias{GenCont_phenotypes}
\docType{data}
\title{
Simulated genetic data with continuous outcomes
}
\description{
Simulated genetic data with continuous outcomes.
}
\usage{data(GenCont)}
\format{
  GenCont is a saved R matrix with 500 rows and 13 columns. The first
  column is the pheotypes and columns 2-13 contain the genotypes. Each
  row represents an indiviaul.
  The same data are stored in flat text files in GenCont_genotypes and GenCont_phenotypes
  (in the directory extdata (in the installed package) or inst/extdata (in the source)).
}

\details{
Genotypes were simulated using FREGENE.
}
\references{
Fregene: Simulation of realistic sequence-level data in populations and ascertained samples Chadeau-Hyam, M. et al, 
2008, BMC Bioinformatics 9:364
}
\examples{
data(GenCont)
}
\keyword{datasets}
