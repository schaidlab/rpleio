pleio.glm.sequential <- function(obj.pleio.glm.fit, pval.threshold){
  pval <- pval.threshold / 2
  n.traits <- obj.pleio.glm.fit$n.traits

  if(all(is.na(obj.pleio.glm.fit$theta))) {
     return(list(pval=NA, count=NA, index.nonzero.coef=NA))
  }
  
  count <- 0
  save <- NULL
  while(pval < pval.threshold & count < n.traits){
    save <- pleio.glm.test(obj.pleio.glm.fit, count.nonzero.coef = count)
    pval <- save$pval
    index.nonzero.coef <- save$index.nonzero.coef
    count <- count + 1
  }

  ## if all traits significant, test is invalid,
  ## return with all traits
  if(count == n.traits & pval <= pval.threshold) {
    index.nonzero.coef <- 1:n.traits
    pval=1.0
  } else {
    ## decrement count to account for "+1" in above loop, in case
    ## pval > pval.threshold when count === 0
    count <- count - 1
  }
  return(list(pval=pval, count=count, index.nonzero.coef=index.nonzero.coef))
}
