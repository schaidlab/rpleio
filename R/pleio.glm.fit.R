pleio.glm.fit  <- function(y, g, glm.family, x.all=NULL, x.index.list=NULL) {

  ## y is a matrix with rows for subjects and cols for traits, with at least 2 traits.
  ## g is a vector of codes for genotype for a SNP, such as 0,1,2, for number of
  ##   minor alleles
  ## glm.family is a vector of character strings for the type of traits in the
  ##     cols of y ("binomial", "gaussian", "ordinal").
  ## x.all is a matrix of covariates for all traits
  ## x.index.list is a vector of lists, where the ith vector item
  ##    (x.index.list[[i]]) is a vector of indices for which cols of x.all
  ##    are used as adjusting covariates for the ith trait in y.

  ## convert y to matrix class
  y <- as.matrix(y)
  
  ## use eps to avoid numerical problems with fitted values
  eps <- 1e-4

  if(ncol(y) < 2){
    stop("less than 2 cols in y")
  }
  
  if( !all(glm.family %in% c("binomial", "gaussian", "ordinal"))  ) {
    stop("unknown glm.family")
  }
  
  n.traits <- ncol(y)
  n.row.y  <- nrow(y)
  
  if(length(glm.family) != n.traits){
    stop("length(glm.family) != n.traits")
  }
  
  if(length(g) != n.row.y){
    stop("length(g) != n.row.y")
  }

  if(is.null(x.all) & !is.null(x.index.list)) {    
    stop("x.all  missing for x.index.list")
  }
  
  ## make a place holder for x.index.list if x.all is NULL
  if(is.null(x.all)){
    x.index.list <- vector(mode="list", length=n.traits)
    for(j in 1:n.traits){
      x.index.list[[j]] <- NA
    }
  }

  if(!is.null(x.all)){
    x.all <- as.matrix(x.all)
    if(is.null(x.index.list)){
      stop("list of X indices not definded")
    }
    if(length(x.index.list) != n.traits){
      stop("length of x.index.list != no. traits")
    }
    
    for(j in 1:n.traits){
      if(any(is.na(x.index.list[[j]]))){
        stop("missing values for x.index.list not allowed")
      }
      if(min(x.index.list[[j]]) < 0){
        stop("x.index.list item < 0")
      }
      if(max(x.index.list[[j]] > ncol(x.all))){
        stop("x.index.list item > ncol(x.all)")
      }
    }
         
    if(nrow(x.all) != nrow(y)){
      stop("number of rows of y and x do not match")
    }
  }


  ## subset to non-missing data for y, X, and  g
  miss.y <- apply(is.na(y), 1, any)

  miss.g <- is.na(g)
    
  no.miss <- !miss.y & !miss.g
  
   if(!is.null(x.all)){
    miss.x <- apply(is.na(x.all), 1, any)
    no.miss <- no.miss & !miss.x
  }
 
  y <- y[no.miss, ]
  g <- g[no.miss]
  if(!is.null(x.all)){
    x.all <- x.all[no.miss,]
  }

  ## if no variance in geno (g), method will fail
  if(length(unique(g))<=1) {
    warning("genotypes (g) are all the same")
    obj <- list(theta=NA, n.intercepts=NA,
              n.coef.covar=NA, n.parm=NA,
              n.traits=ncol(y), an.mat=NA)  
    class(obj) < "pleio.glm.fit"
    return(obj) 
  }

  n.subjs <- nrow(y)
            
  ## check if ordinal traits are numeric, and if so, recode to
  ## have sequential coding 1, 2, ...
  
  for(j in 1:n.traits){
    if(glm.family[j] == "ordinal")
      {
        if(!is.numeric(y[,j])){
          stop( paste("ordinal trait[", j, "] is not numeric", sep="") )
        }     
        y[,j] <- as.numeric(factor(y[,j]))
      }
  }

  n.levels     <- numeric(n.traits)
  n.intercepts <- numeric(n.traits)

  n.traits.expanded <- 0
  
  for(j in 1:n.traits){
    
    if(glm.family[j] == "ordinal"){
      ## for ordinal traits with K levels, there are K-1
      ## expanded traits for the K-1 cumsums
      n.levels[j] <- length(unique(y[,j])) - 1
      n.traits.expanded <- n.traits.expanded  + n.levels[j]
      n.intercepts[j] <- n.levels[j]
    } else {
      n.levels[j] <- 1
      n.traits.expanded <- n.traits.expanded + 1
      n.intercepts[j] <-  1
    }
  }

  intercept <- NULL
  gamma <- numeric(n.traits)
  trait.type.index <- NULL
  
  trait.mat <- NULL
  fit.mat <- NULL
  
  deriv.intercept.mat <- NULL
  deriv.beta.mat <- NULL
  
  fit.stats <- vector(mode="list", length=n.traits)
  
  for(j in 1:n.traits){
    
    if(glm.family[j] =="binomial"){
      if(!any(is.na(x.index.list[[j]])) && any(x.index.list[[j]] %in% 1:ncol(x.all))) {
        x.cov <- as.matrix(x.all[, x.index.list[[j]]])
        fit.stats[[j]]$fit <- glm(y[,j] ~ x.cov  + g, family=stats::binomial)
      } else {
        fit.stats[[j]]$fit <- glm(y[,j] ~  g, family=stats::binomial)
      }
 
      
      fitted.values <- fit.stats[[j]]$fit$fitted.values
      ## patch: avoid fitted = 0 or 1
      fitted.values <- ifelse(fitted.values < eps, eps, fitted.values)
      fitted.values <- ifelse(fitted.values > (1-eps), (1-eps), fitted.values)
      
      trait.mat <- cbind(trait.mat, y[,j])
      fit.mat <- cbind(fit.mat, fitted.values)
      fit.stats[[j]]$var.func <- fitted.values*(1-fitted.values)
      
    } else if(glm.family[j] == "gaussian"){

      if(!any(is.na(x.index.list[[j]])) && any(x.index.list[[j]] %in% 1:ncol(x.all))){
        x.cov <- as.matrix(x.all[, x.index.list[[j]]])   
        fit.stats[[j]]$fit <- glm(y[,j] ~ x.cov  + g, family=stats::gaussian)
      } else {
        fit.stats[[j]]$fit <- glm(y[,j] ~  g, family=stats::gaussian)
      }
          
      fitted.values <- fit.stats[[j]]$fit$fitted.values
      trait.mat <- cbind(trait.mat, y[,j])
      fit.mat <- cbind(fit.mat, fitted.values)
      fit.stats[[j]]$var.func <- rep(1, length(fitted.values))
      
    } else if(glm.family[j] == "ordinal") {

      if(!any(is.na(x.index.list[[j]])) && any(x.index.list[[j]] %in% 1:ncol(x.all))){
        x.cov <- as.matrix(x.all[, x.index.list[[j]]])
        fit.stats[[j]]$fit <- lrm(y[,j] ~ x.cov  + g)
      } else {
        fit.stats[[j]]$fit <- lrm(y[,j] ~  g)
      }
      ## assume predict is predict.glm from stats,
      ## there is a Predict function in rms for lrm that does not give the same options
      ## predict.glm results a matrix
      fitted.values <-  stats::predict(fit.stats[[j]]$fit, type="fitted")
      ## patch: avoid fitted = 0 or 1
      fitted.values <- ifelse(fitted.values < eps, eps, fitted.values)
      fitted.values <- ifelse(fitted.values > (1-eps), (1-eps), fitted.values)
      fit.stats[[j]]$fitted.values <- fitted.values

      ## for ordinal trait, convert from y=j (j=1...K) to
      ## cumsum of indicators of level, excluding 1st cumsum
      y.ord <- y[,j]
      n.unique.levels <- length(unique(y.ord))
      ytemp.mat <- 1*outer(y.ord, 1:n.unique.levels, ">=")
      ytemp.mat <- ytemp.mat[,-1]
      trait.mat <- cbind(trait.mat, ytemp.mat)
      fit.mat <- cbind(fit.mat, fitted.values)
      
      fit.stats[[j]]$var.func <-  fitted.values * (1-fitted.values)
  
    }

  }

  coef.intercept <- NULL
  coef.beta <- NULL
  coef.gamma <- NULL
  n.beta <- numeric(n.traits)
 
  for(j in 1:n.traits){
    coef.temp <- fit.stats[[j]]$fit$coef
    coef.intercept <- c(coef.intercept, coef.temp[1:n.intercepts[j]])
    n.coef <- length(coef.temp)
    n.beta[j] <- n.coef - n.intercepts[j] - 1
    if(n.beta[j] > 0){
      coef.beta <- c(coef.beta, coef.temp[(n.intercepts[j]+1):(n.intercepts[j]+n.beta[j])])
    }
    coef.gamma <- c(coef.gamma, coef.temp[n.intercepts[j] + n.beta[j] + 1])
  }
  
  tot.beta <- sum(n.beta)
  if(tot.beta > 0){
    stop.beta <- cumsum(n.beta)
    start.beta <- stop.beta[1:(length(n.beta)-1)]
    start.beta <- c(1, start.beta + 1)
  }
  
  res.mat <-  trait.mat - fit.mat
  var.func <- NULL

  for(j in 1:n.traits){
    var.func <- cbind(var.func, fit.stats[[j]]$var.func)
  }

  ## matrix of residuals
  res.mat <- res.mat / sqrt(var.func)
  phi <- apply(res.mat^2, 2, mean)
  names(phi) <- NULL
  
  ## correlation matrix
  rmat <- var(t( t(res.mat)/sqrt(phi) )) *(n.subjs-1)/n.subjs

  n.parm <-  length(coef.intercept) + length(coef.beta) + length(coef.gamma)

 
  col.start <- col.end <- numeric(length(n.levels))   
  for(j in 1:length(n.levels))
    {
      if(j == 1){
        col.start[j] <- 1
        col.end[j] <- n.levels[j]
      } else {
        col.start[j] <- col.end[j-1] + 1
        col.end[j] <- col.start[j] + n.levels[j] - 1
      }
    }
  
  start.inter <- stop.inter <- numeric(length(n.levels))
  n.intercepts <- n.levels
  for(j in 1:length(n.intercepts))
    {
      if(j == 1){
        start.inter[j] <- 1
        stop.inter[j] <- n.intercepts[j]
      } else {
        start.inter[j] <- stop.inter[j-1] + 1
        stop.inter[j] <- start.inter[j] + n.intercepts[j] - 1
      }
    }

  
  ## compute an.mat
  
  an.mat <- matrix(0, n.parm, n.parm)
  tot.intercepts <- sum(n.intercepts)
  tot.beta <- sum(n.beta)
 
  for(i in 1:n.subjs){

     deriv.intercept.mat <- matrix(0, tot.intercepts, n.traits.expanded)
     if(tot.beta > 0){
       deriv.beta.mat <-      matrix(0, tot.beta,         n.traits.expanded)
     }
     deriv.gamma.mat <-     matrix(0, n.traits,       n.traits.expanded)
   
     vi <-  numeric(n.traits.expanded)

     for(j in 1:n.traits){
      
       if(glm.family[j] == "binomial"){
         ## binomial
         fval <- fit.stats[[j]]$fit$fitted.values[i]
         deriv.intercept <- fval*(1-fval)
         deriv.intercept.mat[start.inter[j], col.start[j]] <- deriv.intercept
         if(n.beta[j] > 0){
           deriv.beta <- deriv.intercept*x.all[i, x.index.list[[j]]]
           deriv.beta.mat[start.beta[j]:stop.beta[j], col.start[j]] <- deriv.beta
         }
         deriv.gamma <- deriv.intercept * g[i]
         deriv.gamma.mat[j, col.start[j]] <- deriv.gamma
         
       } else if(glm.family[j] == "gaussian"){
    
         ## gaussian
         deriv.intercept.mat[start.inter[j], col.start[j]] <- 1
         if(n.beta[j] > 0){
           deriv.beta <- x.all[i, x.index.list[[j]]]
           deriv.beta.mat[start.beta[j]:stop.beta[j], col.start[j]] <- deriv.beta
         }
         deriv.gamma <-  g[i]
         deriv.gamma.mat[j, col.start[j]] <- deriv.gamma
       } else if(glm.family[j] == "ordinal") {
         ## ordinal
         fval <- fit.stats[[j]]$fitted.values[i,]
         deriv.intercept <- fval*(1-fval)
         deriv.intercept.mat[start.inter[j]:stop.inter[j], col.start[j]:col.end[j]] <- diag(deriv.intercept)
         if(n.beta[j] > 0){
           deriv.beta <- x.all[i, x.index.list[[j]]] %o% deriv.intercept
           deriv.beta.mat[start.beta[j]:stop.beta[j], col.start[j]:col.end[j]] <- deriv.beta
         }
         deriv.gamma <- deriv.intercept * g[i]
         deriv.gamma.mat[j, col.start[j]:col.end[j]] <- deriv.gamma
       }
    
       deriv.mat <- rbind(deriv.intercept.mat, deriv.beta.mat, deriv.gamma.mat)
       if(glm.family[[j]] == "ordinal"){
         vi[col.start[j]:col.end[j]] <- fit.stats[[j]]$var.func[i,]
       } else {
         vi[col.start[j]:col.end[j]] <- fit.stats[[j]]$var.func[i]
       }

     }

     se <- diag(sqrt(phi * vi))
     vari <- se %*% rmat %*% se
     vinv <- solve(vari)
     
     ## an is a matrix of deriv * vinv * deriv, where
     ## deriv is deriv mu wrt theta
     
     an.mat <- an.mat + deriv.mat %*% vinv %*% t(deriv.mat)
   }
     
  theta <- NULL
 
  ## collect parm estimates
  theta <- coef.intercept
  if(tot.beta > 0) {
    theta <- c(theta, coef.beta)
  }
  theta <- c(theta, coef.gamma)
    
  obj <- list(theta=theta,
              n.intercepts=tot.intercepts,
              n.coef.covar=tot.beta,
              n.parm=n.parm,
              n.traits=n.traits,
              an.mat=an.mat)
  
  class(obj) < "pleio.glm.fit"

  return(obj)
       
}
