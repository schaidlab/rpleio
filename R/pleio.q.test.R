## Pleiotropy for quantitative trait, test Ho of count.nonzero.beta.
## To reject Ho is to say at least one more are significant.
## Author: DJ Schaid
## Date: 7/27/2016

pleio.q.test <- function(obj.fit, count.nonzero.beta = 0){
  if(all(is.na(obj.fit$x))) {
     return(list(stat=NA, pval=NA, df=NA, index.nonzero.beta=NA, tests=NA))
  }
  x <- obj.fit$x
  xx.inv <- obj.fit$xx.inv
  beta.ols <- obj.fit$beta.ols
  n.traits <- obj.fit$n.traits
  
  if(count.nonzero.beta < 0){
    stop("count.nonzero.beta < 0")
  }
  if(count.nonzero.beta > (n.traits - 1)){
    stop("count.nonzero.beta > (n.traits - 1)")
  }
  
  if(count.nonzero.beta == 0)
    {
      vk <- diag(n.traits)
      beta.star <- xx.inv %*% t(vk) %*% solve(vk %*% xx.inv %*% t(vk)) %*% vk %*% beta.ols
      temp <- x %*% beta.star
      stat <- sum(temp^2)
      tk <- stat
      lrt <- stat
      df <- n.traits
      pval  <- 1-pchisq(lrt, df=df)
      index.nonzero.beta <- 0
      vk.set <- 0
    }
  
  if(count.nonzero.beta > 0)
    {
      
      vk.set <- combn(n.traits, count.nonzero.beta)
      
      ## rows of vk.set are the set of indices for beta's that
      ## are not constrained to 0, so these are the indices
      ## used to create vk for testing vk*beta = 0.
      ## The cols of vk.set are different possible configurations,
      ## requiring computation of tk for each col of indices, then
      ## the final lrt = min(tk)
      
      set.size <- ncol(vk.set)
      tk <- numeric(set.size)
      index.constrained <- 1
      tmin <- 10^6
      
      for(j in 1:set.size){
        
        vk <- diag(n.traits)
        index.exclude <- vk.set[,j]
        vk <- vk[-index.exclude,, drop=FALSE]
        beta.star <- xx.inv %*% t(vk) %*% solve(vk %*% xx.inv %*% t(vk)) %*% vk %*% beta.ols
        temp <- x %*% beta.star
        
        tk[j] <- sum(temp^2)
        
        if(tk[j] < tmin){
          tmin <- tk[j]
          index.constrained <- j
        }
      }
      
      lrt <- min(tk)
      df <- n.traits - nrow(vk.set)
      pval <- 1 - pchisq(lrt, df=df)
      index.nonzero.beta <- vk.set[, index.constrained]
    }


  free.betas <- cbind(t(vk.set), tk)
  free.betas <- data.frame(free.betas)
  n.freebeta <- nrow(vk.set)

  ## patch for when count.nonzero.beta == 0
  if(is.null(n.freebeta)) n.freebeta <- 1
 
  names(free.betas)<- c(paste("index", 1:n.freebeta, sep="."), "tk")

  return(list(stat=lrt, pval=pval, df=df, index.nonzero.beta=index.nonzero.beta,
              tests=free.betas))
   
}
