
## Pleiotropy for quantitative trait, sequential test of 0, then 1, etc.
## Author: DJ Schaid
## Date: 7/27/2016
pleio.q.sequential <- function(obj.fit, pval.threshold){
  pval <- pval.threshold / 2
  n.traits <- obj.fit$n.traits

  if(all(is.na(obj.fit$x))) {
    return(list(pval=NA, index.beta=NA))
  }
  
  count <- 0
  save <- NULL
  while(pval < pval.threshold & count < n.traits){
    save <- pleio.q.test(obj.fit, count.nonzero.beta=count)
    pval <- save$pval
    index.beta <- save$index.nonzero.beta
    count <- count + 1
  }

  ## if all traits significant, test is invalid,
  ## return with all traits
  if(count == n.traits & pval <= pval.threshold) {
    index.beta <- 1:n.traits
    pval=1.0
  } else {
    ## last trait added not signif, so 
    ## decrement count to account for "+1" in above loop, in case
    ## pval > pval.threshold when count === 0
    count <- count - 1
  }
    
  return(list(pval=pval, index.beta=index.beta))
}
