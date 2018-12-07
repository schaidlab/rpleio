\name{pleio.qdemo}
\alias{pleio.qdemo}
\docType{data}
\title{
 Demonstration dataset for quantitative pleiotropy tests
}
\description{
A simulated dataset with 6 y variables simulated from multivariate
normal distribution with common correlation structure, correlation of
0.5, and genotypes simulated based on minor allele frequency of 0.2, and
assume that betas for traits 2 and 3 have non-zero values, while all
other traits are not asssociated with dose of minor allele. 
}
\usage{data("pleio.qdemo")}
\format{
  y is a data.frame with 500 observations (rows) and 6 quantitative traits (columns),
  while geno is a single SNP of dosage of a minor allele,
  simulated from frequency of 0.2. 
}
\examples{
data(pleio.qdemo)
str(y)
table(geno)
}
\keyword{datasets}
