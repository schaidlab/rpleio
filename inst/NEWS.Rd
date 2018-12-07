\name{NEWS}
\title{NEWS file for the pleio package}
\section{Changes in version 1.5}{
  \itemize{
    \item Add new functions for GLMs to extend the same functionality of
    pleio.q.fit, pleio.q.test, and pleio.q.sequential to functions of
    the same name, with the q replaced with glm.  The pleio.glm.fit
    function requires the families of the traits (y) and indices of
    covariates from a user-defined x matrix, and vectors of indices of
    the x matrix that are to be modeld with each trait (y).
    \item added vignette for pleio.glm.*, with similar steps as pleio.q.*
    \item If all traits are significant, return pval=1 and indices of
    all traits.
    \item Allow both matrix and data.frame objects for covariates (x.all)
    \item If all genotypes are the same, return empty objects.
    \item Allow either matrix or data.frame for traits (fix for v1.6)
  }
}
\section{Changes in version 1.3}{
  \itemize{
    \item Initial CRAN release
    \item Include main functions pleio.q.fit, pleio.q.test and
    pleio.q.sequential, as explained in the vignette pleioQ
  }
}
