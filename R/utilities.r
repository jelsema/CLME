
##
## Various support functions to assist the main functions
##

# Function to define the model matrices for clme
model_terms_clme <- function( formula, data, ncon=1 ){
  
  formula2 <- update.formula( formula , . ~ . - 1 )
  
  nn    <- nrow(data)
  Unull <- c( rep("RemoveAAA",round(nn/2)), rep("RemoveBBB",nn-round(nn/2)) )
  
  data  <- cbind( data , Unull )
  
  formula3 <- update.formula( formula , . ~ . + (1|Unull) - 1)
  
  suppressMessages( mterms   <- lme4::lFormula( formula3, data=data ) )
  
  Y  <- mterms$fr[,1]
  X  <- mterms$X
  P1 <- sum( attr(X, "assign") <= ncon )
  X1 <- X[,     1:P1     , drop=FALSE]
  
  if( ncol(X) > P1 ){
    X2 <- X[,(P1+1):ncol(X), drop=FALSE]  
  } else{
    X2 <- NULL
  }
  
  U            <- t( as.matrix(mterms$reTrms$Zt) )
  dframe       <- mterms$fr
  dframe$Unull <- NULL
  
  if( ncol(U)==2 ){
    mmat <- list( Y=Y, X1=X1, X2=X2, P1=P1, U=NULL, formula=formula2, dframe=dframe, 
                  REidx=NULL, REnames=NULL )
  } else{
    drop.col <- which( colnames(U) %in% c("RemoveAAA","RemoveBBB") )
    drop.nam <- which( colnames(mterms$reTrms$flist)=="Unull" )
    
    U <- U[, -drop.col]
    
    REnames <- colnames(mterms$reTrms$flist)[-drop.nam]
    REidx   <- mterms$reTrms$Lind[-drop.col]
    
    mmat <- list( Y=Y, X1=X1, X2=X2, P1=P1, U=U, formula=formula2, dframe=dframe, 
                  REidx=REidx, REnames=REnames )
  }
  
  return( mmat )
  
}


##
## Some methods for class CLME
##



is.clme <- function(x) inherits(x, "clme")




as.clme <- function( x , ... ){
  
  if( is.clme(x) ){
    return(x)
  } else{
    
    err.flag  <- 0
    flagTheta <- flagSsq <- flagTsq <- flagCov <- flagW1 <- flagW2 <- flagP1 <- flagP2 <- flagConst <- ""
    
    if( !is.numeric(x$theta) ){
      err.flag  <- 1
      flagTheta <- " theta must be numeric \n"
      x$theta   <- numeric(0)
    }
    
    if( !is.numeric(x$ssq) ){
      err.flag <- 1
      flagSsq  <- " ssq must be numeric \n"
      x$ssq    <- numeric(0)
    } 
    
    if( !is.null(x$tsq) & !is.numeric(x$tsq) ){
      err.flag <- 1
      flagTsq  <- " if present, tau must be numeric \n"
      x$tsq    <- NULL
    }
    
    if( !is.matrix(x$cov.theta) || !is.numeric(x$cov.theta) ||
          nrow(x$cov.theta) != ncol(x$cov.theta) ||
          nrow(x$cov.theta) != length(x$theta)   ||
          sum(sum(abs(x$cov.theta - t(x$cov.theta)))) > sqrt(.Machine$double.eps) ){
      err.flag    <- 1
      flagCov     <- " cov.theta must be square, symmetric, numeric matrix with dimensions equal to length of theta\n"
      x$cov.theta <- matrix( numeric(0) , nrow=length(x$theta) , ncol=length(x$theta) )
    }
    
    if( !is.numeric(x$ts.glb) ){
      err.flag <- 1
      flagW1   <- " ts.glb must be numeric \n"
      x$ts.glb  <- numeric(0)
    } 
    
    if( !is.numeric(x$ts.ind) ){
      err.flag <- 1
      flagW2   <- " ts.ind must be numeric \n"
      x$ts.ind  <- numeric(0)
    }
    
    if( !is.numeric(x$p.value) || length(x$p.value) != length(x$ts.glb) ){
      err.flag   <- 1
      flagP1     <- " p.value must be numeric and of same length as ts.glb \n"
      x$p.value  <- numeric(0)
    } 
    
    if( !is.numeric(x$p.value.ind) || length(x$p.value.ind) != length(x$ts.ind) ){
      err.flag       <- 1
      flagP2         <- " p.value.ind must be numeric and of same length as ts.ind \n"
      x$p.value.ind  <- numeric(0)
    } 
    
    if( !is.list(x$constraints) ){
      err.flag        <- 1
      flagConst       <- " constraints must be list \n"
      x$constraints   <- list( A=matrix( numeric(0) ) )
    } else{
      cnames <- names(x$constraints)
      if( sum(cnames=="A") != 1 ){
        err.flag        <- 1
        flagConst       <- " constraints must contain element A\n"
        x$constraints$A <- matrix( numeric(0) , nrow=length(x$ts.ind ))
      }
    }
    
    
    if( err.flag==1 ){
      err.mssg <- paste( "coercing 'x' to class 'clme' produced errors: \n", 
                         flagTheta, flagSsq, flagTsq, flagCov, flagW1,
                         flagW2, flagP1, flagP2, flagConst, "output may not be valid." , sep = "")
      # warning(warn, sys.call(-1))
      warning( err.mssg )
    }
    
    class(x) <- "clme"
    return(x)
    
  }
}


################################################################################


AIC.clme <- function( object, ..., k=2 ){
  ## For BIC, set k = ln( n/(2*pi) )
  logl <- logLik.clme( object, ...)[1]
  kk   <- ncol(model.matrix.clme(object)) + length(object$tsq) + length(object$ssq) 
  # aic  <- 2*(kk - logl)
  aic  <- k*kk - 2*logl
  return(aic)
}


confint.clme <- function(object, parm, level=0.95, ...){
  ## More types of confidence intervals (e.g., bootstrap) may be added in the future.
  ## If so, default confidence interval will be the current methods
  
  ## Confidence intervals based on unconstrained variance-covariance matrix
  ## Actual covarage >= Nominal coverage
  
  cc     <- match.call()
  digits <- cc$digits
  if( is.null(digits) ){     digits <- 3 }
  if( !is.numeric(digits) ){ digits <- 3 }
  if( digits < 0  ){         digits <- 3 }
  
  
  alpha <- 1 - level
  theta <- fixef(object)
  cv    <- qnorm(1-alpha/2)
  varco <- vcov( object )
  lcl   <- as.numeric( format( round(theta - cv*sqrt(diag(varco)) , digits=digits)) )
  ucl   <- as.numeric( format( round(theta + cv*sqrt(diag(varco)) , digits=digits)) )
  
  ints <- cbind( lcl, ucl)
  
  ## Return intervals
  return( ints )
}


fixef.clme <- function( object, ... ){
  ## Print out the fixed effects
  if( is.clme(object) ){
    return( object$theta )
  } else{
    stop("'object' is not of class clme")
  }
}
fixed.effects.clme <- function( object, ... ){
  fixef.clme( object, ... )
}
#coefficients.clme <- function( object, ... ){
#  fixef.clme( object, ... )
#}
#coef.clme <- function( object, ... ){
#  fixef.clme( object, ... )
#}


formula.clme <- function(x, ...){
  return( x$formula )
}


logLik.clme <- function( object, ...){
  ## Residuals
  YY <- object$dframe[,1]
  XX <- model.matrix( object )
  TT <- fixef.clme( object )
  RR <- YY - apply( XX , 1 , FUN=function(xx,tht){ sum(xx*tht) }, tht=TT )
  nn <- nobs.clme(object)
  
  ## Covariance matrix (piecewise)
  ssq <- object$ssq
  Nks <- object$gfix
  ssqvec <- rep( ssq, Nks )
  RSiR <- c( t(RR) %*% (RR/ssqvec) )
  detS <- sum( log(ssqvec) )
  
  if( is.null(object$Qs) ){
    ## Fixed effects only
    RPiR   <- RSiR
    detPhi <- detS
  } else{
    ## Mixed Effects    
    mmat   <- model_terms_clme( object$formula, object$dframe )
    UU     <- mmat$U
    tsq    <- object$tsq
    Qs     <- object$gran
    tsqvec <- rep( tsq, Qs )
    RSiU   <- matrix( apply( UU, 2, FUN=function(uu,rr){ sum(uu*rr) }, rr=(RR/ssqvec) ), nrow=1 )
    
    U1         <- apply( UU, 2, FUN=function(uu,sq){uu/sq}, sq=sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve( tusu )
    
    RPiR   <- RSiR - c( RSiU%*%(tusui%*%t(RSiU)) )
    detPhi <- det( tusui ) * sum( log(tsqvec) )  * detS
  }
  
  logL <- c( (-nn*length(TT)*log(2*pi)) - detPhi - RPiR )/2
  return( logL )
  
}



model.frame.clme <- function( formula , ...){
  ## Return the data frame
  mmat     <- model_terms_clme( formula, ... )  
  return( mmat$dframe )
}

model.matrix.clme <- function( object, ...){
  ## Return the fixed-effects matrix
  mmat <- model_terms_clme( object$formula, object$dframe )
  X1   <- mmat$X1
  X2   <- mmat$X2
  return( cbind(X1, X2) )  
}

nobs.clme <- function(object, ...){
  nrow( model.matrix.clme(object) )
}

print.clme <- function(x, ...){
  #cc     <- match.call()
  #digits <- cc$digits
  
  object <- x
  ## Print out residuals of specified type
  if( !is.clme(object) ){
    stop("'object' is not of class clme")
  }
  
  cat( "Linear mixed model subject to order restrictions\n")
  cat( "Formula: ")
  print( object$formula )  
  crit <- c(logLik.clme(object),
            AIC.clme(object),
            AIC( object, k=log(nobs.clme(object)/(2*pi)) ) )  
  
  critc <- format( crit , digits=5)
  
  cat( "\nlog-likelihood:", critc[1] )
  cat( "\nAIC:           ", critc[2] )
  cat( "\nBIC:           ", critc[3] )
  cat( "\n(log-likelihood, AIC, BIC computed under normality)")  
  
  cat( "\n\nFixed effect coefficients (theta): \n")
  print( fixef.clme(object)  )
  
  cat( "\nVariance components: \n")
  print( VarCorr.clme(object) )
  
  cat( "\n\nModel based on", object$nsim, "bootstrap samples." )
  
}


ranef.clme <- function( object, ... ){
  ## Print out the random effects
  if( is.clme(object) ){
    return( object$random.effects )
  } else{
    stop("'object' is not of class clme")
  }
}
random.effects.clme <- function( object , ... ){
  ranef.clme( object, ... )
}


residuals.clme <- function( object, type="FM", ... ){
  ## Print out residuals of specified type
  if( is.clme(object) ){
    ridx <- which( c("PA", "SS", "FM")==type )
    if( ncol(object$residuals)<ridx ) ridx <- 1
    return( object$residuals[,ridx] )
  } else{
    stop("'object' is not of class clme")
  }
}


sigma.clme <- function( object, ...){
  return( object$ssq )
}


VarCorr <- function( object, ...){ UseMethod("VarCorr") }

VarCorr.clme <- function(object, ...){
  ## Print out variances or SDs
  ## Defines tiny class "varcorr_clme" to handle printing
  ## using the method: print.varcorr_clme
  if( !is.clme(object) ){
    stop("'object' is not of class clme")
  } else{
    varcomps <- matrix( c(object$tsq, object$ssq ), ncol=1 )
    rnames   <- c( "Source", names(object$tsq), names(object$ssq) )
    rownames(varcomps) <- rnames[-1]
    colnames(varcomps) <- "Variance"
    class(varcomps)    <- "varcorr_clme"
    return( varcomps )
  } 
}

## Leave this method out of alphabetical order so that
## is it right next to the VarCorr.clme method
print.varcorr_clme <- function( x, ... ){
  varcomps <- x
  rnames   <- c( "Source", rownames( varcomps ) )
  rnames   <- str_pad(rnames, width=max(nchar(rnames)), side = "right", pad = " ")
  vars     <- format( varcomps , digits=5)
  
  cat( rnames[1], "\t" , "Variance" )
  for( ii in 1:length(vars) ){
    cat( "\n", rnames[ii+1], "\t" , vars[ii] )
  }
}


vcov.clme <- function(object, ...){
  ## Print out covariance matrix of theta
  if( is.clme(object) ){
    return( object$cov.theta )
  } else{
    stop("'object' is not of class clme")
  }  
}


##
## Hidden functions to format characters / decimals
##



## Align the length of header with length of values
.align_table.clme <- function( tble, digits=4, ... ){
  cnames <- colnames( tble )
  for( ii in 1:length(cnames) ){
    ntitle <- nchar( cnames[ii] )
    maxc   <- max( nchar(tble[,ii]) )    
    if( ntitle > maxc ){
      tble[,ii] <- str_pad( tble[,ii], width=ntitle, side = "left", pad = " ")    
    }
    if( ntitle < maxc ){
      cnames[ii] <- str_pad( cnames[ii], width=(maxc+1), side = "right", pad = " ")    
    }
  }
  colnames(tble) <- cnames
  return( tble )
}





