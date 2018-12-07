\name{pleio.q.fit}
\alias{pleio.q.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Intermediate Statisics for Testing Pleiotropy
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Pre-compute statisics for testing pleiotropy, to be used in
  pleio.q.test and pleio.q.sequential. This is an earlier function for
  gaussian traits without covariates. The function pleio.glm.fit is more
  general, allowing for gaussian traits, and possibly adjusting covariates. 
}
\usage{
pleio.q.fit(y, geno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ matrix of multiple quantitative traits to test for
    pleiotropy on genotype vector. Rows are for subjects and columns are for traits.
}
  \item{geno}{vector of genotypes coded as dosage of a minor allele,
    0,1,2. Allows imputed genotype dosage.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The matrix y of traits is expected to have  rows for subjects and
  columns for traits. The number of rows of y  should match the length of
  vector geno.  This function computes intermediate statistics that are
  subsequently used to test various hypotheses about pleiotropy.
  For technical details,  Schaid et al. (Genetics, 2016)
}
\value{
%%  ~Describe the value returned
A list containing
  \item{x}{the de-correlated genotype design matrix}
  \item{xx.inv}{inverse of x'x}
  \item{beta.ols}{estimated unconstrained beta coefficients for y traits based on
    ordinary least squares}
  \item{n.traits}{number of traits (columns) in y}
}
\references{
  %% ~put references to the literature/web site here ~
  Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. Genetics. 2016
  Oct;204(2):483-497. 
}
\author{
  Dan Schaid, Jason Sinnwell.
%%  ~~who you are~~
}

