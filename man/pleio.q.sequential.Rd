\name{pleio.q.sequential}
\alias{pleio.q.sequential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform sequential tests of pleiotropy
}
\description{
  Perform sequential tests of pleiotropy to determine the number of
  traits, and which traits, are associatd with a vector of
  genotypes. For gaussian traits without covariates.
}
\usage{
pleio.q.sequential(obj.fit, pval.threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.fit}{ result of pleio.q.fit }
  \item{pval.threshold}{p-value for rejecting the null hypothesis of the
    specified number of non-zero betas.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Perform sequential tests of pleiotropy, starting at the usual
  multivarite null hypothesis that all coefficients = 0. If this test
  rejects because the p-value < pval.threshold, then allow one coefficient to
  be non-zero in order to test whether the remaining coefficients = 0.
  If the test of one non-zero coefficient rejects, then allow two
  non-zero coefficients, considering all possible combinations of two 
  non-zero coefficients and test whether the remaining coefficients =
  0. Continue this sequential testing until the p-value for a test is
  greater than the specific pval.threshold. The step at which the
  p-value > pval.threshold determines which traits are associated with
  the genotype. If there are m traits, the sequential testing stops
  either when p-value > pval.threshold, or when (m-1) traits are tested.
  If the p-value remains less than pval.threshold when testing (m-1)
  traits, this implies that all m traits are associated with the
  genotype.
}
\value{
%%  ~Describe the value returned
A list containing:
\item{pval }{p-value of the final test from the sequential testing}
\item{index.beta}{index of columns of y that have non-zero betas.
  These indices indicate which traits are associated with the genotype,
  accounting for the correlations among the traits.}
}
\references{
  %% ~put references to the literature/web site here ~
   Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. To appear in Genetics, 2016.
}
\author{
  %%  ~~who you are~~
  Dan Schaid and Jason Sinnwell
}
\examples{
data(pleio.qdemo)
fit <- pleio.q.fit(y, geno)
test.seq <- pleio.q.sequential(fit, pval.threshold=.05)
test.seq
}

