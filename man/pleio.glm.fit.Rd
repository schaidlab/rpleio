\name{pleio.glm.fit}
\alias{pleio.glm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Intermediate Statisics for Testing Pleiotropy
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Pre-compute statisics for testing pleiotropy, to be used in
  pleio.glm.test and pleio.glm.sequential 
}
\usage{
pleio.glm.fit(y, g, glm.family, x.all=NULL, x.index.list=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ matrix of multiple traits to test for
    pleiotropy on genotype vector. Rows are for subjects and columns are for traits.
  }
  \item{g}{vector of genotypes coded as dosage of a minor allele, 0,1,
    and 2. Allows imputed genotype dosage.
  }
  \item{glm.family}{vector of strings corresponding to the type of
    phenotype for the columns of y. Types accepted are gaussian,
    binomial, and ordinal.
  }
  \item{x.all}{ matrix of all covariates traits to adjust for the traits.
  }
  \item{x.index.list}{ a list of vectors containing the index of the
    columns of x.all for which to be adjusted for each trait. The list
    must have the same number of vectors as there are traits in y. If
    one trait has no covarites, specify a vector containing 0 for that trait.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The matrix y of traits is expected to have  rows for subjects and
  columns for traits. The number of rows of y  should match the length of
  vector g. If adjusting traits by covarites in x.all, the x.index.list
  must be a list containing a vector of indices of the columns of x.all
  for each trait. If one trait does not have covariates and others do,
  then the vector for the trait with no covariate must be a vector containing
  the number 0. This function computes intermediate statistics that are
  subsequently used to test various hypotheses about pleiotropy.
  For technical details,  Schaid et al. (2016, 2017)
}
\value{
%%  ~Describe the value returned
A list containing
  \item{theta}{  }
  \item{n.intercepts}{ }
  \item{n.coef.covar}{ }
  \item{n.parm}{ }
  \item{an.mat}{ }
  \item{n.traits}{number of traits (columns) in y}
}
\references{
  %% ~put references to the literature/web site here ~
  Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. Genetics. 2016
  Oct;204(2):483-497. 

  Schaid DJ, Tong X, Batzler A, Sinnwell JP, Qing J, Biernacka
  JM. Multivariate Generalized Linear Model for Genetic
  Pleiotropy. Under review. 
}
\author{
  Dan Schaid, Jason Sinnwell.
%%  ~~who you are~~
}

