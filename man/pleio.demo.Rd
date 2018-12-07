\name{pleio.demo}
\alias{pleio.demo}
\alias{geno}
\alias{y}
\alias{x}
\docType{data}
\title{
 Demonstration dataset for pleiotropy tests with a mixture of trait
 types and covariates
}
\description{
A simulated dataset with 3 y variables simulated from multivariate
normal with common correlation structure, correlation of 0.5. The traits
are (1) gaussian, (2) binomial, and (3) ordinal. Genotypes were simulated based
on minor allele frequency of 0.2, and assume that the 3 traits 
are not asssociated with dose of minor allele (all 3 betas = 0).
}
\usage{data("pleio.demo")}
\format{
  y is a data.frame with 500 observations and 3 traits, whose
  distribution families are gaussian, binomial, and ordinal, respectively;  
  x is a matrix of 5 covariates,
  geno is a single SNP of dosage of a minor allele simulated from frequency of 0.2. 
}
\examples{
data(pleio.demo)
str(y)
table(geno)
x[1:5,]
}
\keyword{datasets}
