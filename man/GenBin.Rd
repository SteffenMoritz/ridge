\name{GenBin}
\alias{GenBin}
\alias{GenBin_genotypes}
\alias{GenBin_phenotypes}
\docType{data}
\title{
Simulated genetic data with a binary phenotypes
}
\description{
Simulated genetic data at 15 SNPs, together with simulated binary phenotypes
}
\usage{data(GenBin)}
\format{
  GenBin is a saved R matrix with 500 rows and 15 columns. The first
  column is the pheotypes and columns 2-15 contain the genotypes. Each
  row represents an indiviaul.
  The same data are stored in flat text files in GenBin_genotypes and GenBin_phenotypes
  (in the directory extdata (in the installed package) or inst/extdata (in the source)).
}
\source{
Simulated using FREGENE
}
\references{
Fregene: Simulation of realistic sequence-level data in populations and ascertained samples Chadeau-Hyam, M. et al, 
2008, BMC Bioinformatics 9:364
}
\examples{
data(GenBin)
}
\keyword{datasets}
