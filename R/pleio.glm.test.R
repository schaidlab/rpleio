pleio.glm.test  <- function(obj.pleio.glm.fit, count.nonzero.coef = 0){

  if(all(is.na(obj.pleio.glm.fit$theta))) {
    return(list(stat=NA, pval=NA, df=NA,
              index.nonzero.coef=NA, vk.set=NA, tk=NA))
  }
  if(count.nonzero.coef >= obj.pleio.glm.fit$n.traits){
    stop("count.nonzero.coef >= no. traits")
  }
  if(count.nonzero.coef < 0){
    stop("count.nonzero.coef < 0")
  }
  
  theta <- obj.pleio.glm.fit$theta
  n.intercepts <- obj.pleio.glm.fit$n.intercepts
  n.coef.covar <- obj.pleio.glm.fit$n.coef.covar
  n.parm <- obj.pleio.glm.fit$n.parm
  n.traits <- obj.pleio.glm.fit$n.traits
  an.mat <- obj.pleio.glm.fit$an.mat
  n.other <- n.intercepts + n.coef.covar

  ## use index.other to indentify parms other than coefs for genetic variant
  ## that will be excluded from constrained test (i.e., intercepts and
  ## betas for adjusting covars)
  
  index.other <- seq(from=1, to=n.other)
  an.inv <- solve(an.mat)
 
  if(count.nonzero.coef == 0)
  {
    vk <- diag(n.parm)
    ## remove intercepts & beta's for covars
    vk <- vk[-index.other, ]
    theta.star <- as.vector(an.inv %*% t(vk) %*% solve(vk %*% an.inv %*% t(vk)) %*% vk %*% theta)
    stat <- as.vector(theta.star %*% an.mat %*% theta.star)
    tk <- stat
    df <- n.traits
    pval  <- 1-pchisq(stat, df=df)
    index.nonzero.coef <- 0
    vk.set <- 0
  }
  
  if(count.nonzero.coef > 0){

    vk.set <- combn(n.traits, count.nonzero.coef)
    
    ## rows of vk.set are the set of indices for gamma's that
    ## are not constrained to 0, so these are the indices
    ## used to create vk for testing vk*gamma = 0.
    ## The cols of vk.set are different possible configurations,
    ## requiring computation of tk for each col of indices, then
    ## the final stat = min(tk)
    
    set.size <- ncol(vk.set)
    tk <- numeric(set.size)
    index.constrained <- 1
    tmin <- 10^6

    for(j in 1:set.size){
      
      vk <- diag(n.parm)
      
      index.exclude <- c(index.other, (n.other + vk.set[,j]))
     
      vk <- vk[-index.exclude,, drop=FALSE] 
      
      theta.star <- as.vector(an.inv %*% t(vk) %*% solve(vk %*% an.inv %*% t(vk)) %*% vk %*% theta)
      stat <- as.vector(theta.star %*% an.mat %*% theta.star)
      tk[j] <- stat
      
      if(tk[j] < tmin){
        tmin <- tk[j]
        index.constrained <- j
      }
    }
    
    stat <- min(tk)
    df <- n.traits - nrow(vk.set)
    pval <- 1 - pchisq(stat, df=df)
    index.nonzero.coef <- vk.set[, index.constrained]
    
  }
  
  return(list(stat=stat,
              pval=pval,
              df=df,
              index.nonzero.coef=index.nonzero.coef,
              vk.set=vk.set,
              tk=tk))
}

