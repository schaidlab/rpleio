\name{pleio.q.test}
\alias{pleio.q.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single test of the number of traits associated with genotype
}
\description{
  Perform single test of the number of traits associated with a
  genotype, by allowing a specified number of regression beta's to
  be unconstrained.
}
\usage{
pleio.q.test(obj.fit, count.nonzero.beta = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.fit}{ result of pleio.q.fit }
  \item{count.nonzero.beta}{Number of non-zero betas. A non-zero beta is
    allowed to be unconstrained, while all other beta's are constrained
    to be zero.
    }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  By specifying the number of non-zero beta's, the algorithm evaluates
  all possible ways of choosing unconstrained and constrained beta's, and
  for each configuration a statistic (tk) is computed. This tk statistic
  can be considered a measure of fit of a model. The minimum tk over all
  possible configurations provides a global test of whether one of the
  models fits well.
}
\value{
%%  ~Describe the value returned
A list containing:
  \item{stat }{global test statistict}
  \item{df}{degrees of freedom of the statistic}
  \item{pval}{p-value for the test}
  \item{index.nonzero.beta}{index of the non-zero beta(s) that provide(s)
    the minimum tk goodness of fit statistic - this configuration is
    assumed to have beta's for all other indices equal to zero.}
  \item{tests}{data.frame containing the tests performed.
    For m traits, and k = count.nonzero.beta, there are
    m-choose-k tests considered in the null hypothesis. The data.frame
    provides the indices of the unconstrained betas and the
    corresponding tk test statistic for the configuration.}
}
\references{
  %% ~put references to the literature/web site here ~
  Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. Genetics. 2016
  Oct;204(2):483-497.
}
\author{
  %%  ~~who you are~~
  Dan Schaid and Jason Sinnwell
}
\examples{
data(pleio.qdemo)
fit <- pleio.q.fit(y, geno)
## usual multivariate test of whether all betas = 0
test0 <- pleio.q.test(fit, count.nonzero.beta = 0)
test0
## test whether allowing 2 betas to be non-zero fits data
test2 <- pleio.q.test(fit, count.nonzero.beta = 2)
test2
}

