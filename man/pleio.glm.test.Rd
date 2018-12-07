\name{pleio.glm.test}
\alias{pleio.glm.test}
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
pleio.glm.test(obj.pleio.glm.fit, count.nonzero.coef = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.pleio.glm.fit}{ result of pleio.glm.fit }
  \item{count.nonzero.coef}{Number of non-zero coefficients (betas). A non-zero beta is
    allowed to be unconstrained, while all other beta's are constrained to be zero.
    }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  By specifying the number of non-zero coefficients, the algorithm evaluates
  all possible ways of choosing unconstrained and constrained betas, and
  for each configuration a statistic (tk) is computed. This tk statistic
  can be considered a measure of fit of a model. The minimum tk over all
  possible configurations provides a global test of whether one of the
  models fits well.
}
\value{
%%  ~Describe the value returned
A list containing:
  \item{stat }{global test statistic}
  \item{df}{degrees of freedom of the statistic}
  \item{pval}{p-value for the test}
  \item{index.nonzero.coef}{index of the non-zero coefficients that provide
    the minimum tk goodness of fit statistic - this configuration is
    assumed to have coefficients for all other indices equal to zero.}
  \item{tk}{test testistic for the tests performed on trait combinations
    in vk.set}
  \item{vk.set}{data.frame containing the tests performed.
    For m traits, and k = count.nonzero.coef, there are
    m-choose-k tests considered in the null hypothesis. The data.frame
    provides the indices of the unconstrained coefficients for the
    corresponding tk test statistic for the configuration. Rows are the
    indices for each configuration, and the columns are for the
    different configurations tested.}
}
\references{
  Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. Genetics. 2016
  Oct;204(2):483-497. 
}
\author{
  %%  ~~who you are~~
  Dan Schaid and Jason Sinnwell
}
\examples{
data(pleio.demo)
obj <- pleio.glm.fit(y, geno, glm.family=c("gaussian","binomial","ordinal"))
test1 <- pleio.glm.test(obj, count.nonzero.coef = 0)
test1
test2 <- pleio.glm.test(obj, count.nonzero.coef = 1)
test2
}

