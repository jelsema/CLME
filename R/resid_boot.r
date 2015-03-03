resid_boot <-
function(formula, data, gfix=NULL, eps=NULL, xi=NULL, null.resids=TRUE,
         theta=NULL, ssq=NULL, tsq=NULL, cov.theta=NULL, seed=NULL, 
         nsim=1000, mySolver="LS", ncon=1, ... ){
  
  ##
  ## I should consider a way to condense this so it doesn't need to 
  ## be copied between here, clme() and clme_resids().
  ##
  suppressMessages( mmat <- model_terms_clme( formula, data, ncon ) )
  formula2 <- mmat$formula
  Y  <- mmat$Y
  P1 <- mmat$P1
  X1 <- mmat$X1
  X2 <- mmat$X2
  U  <- mmat$U
  
  if( is.null(gfix) ){
    gfix <- rep("Residual", nrow(U)) 
  }
  Nks <- table(gfix)
  
  if( !is.null(U) ){
    if( !is.null(mmat$REidx) ){
      Qs        <- table( mmat$REidx )
      names(Qs) <- mmat$REnames
    } else{
      Qs  <- table( rep("tsq", ncol(U)) )
    }    
  } else{
    Qs <- NULL
  }
  #############################################
  
  if( is.numeric(seed) ){
    set.seed(seed)
  }
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix( cbind( X1,X2 ))
  K  <- length(Nks)
    
  ## Estimate parameters if they are missing
  
  if( is.null(ssq) ){
    theta1 <- ginv( t(X)%*%X )%*%(t(X)%*%Y)  
    ssq   <- vector()
    for( k in 1:K ){
      Yk <- Y[ N1[k]:N2[k] ]
      Xk <- X[ N1[k]:N2[k],]  
      ssq[k] <- sum( (Yk - Xk%*%theta1)^2 ) / (Nks[k])
    }
  }
  
  if( (Q > 0) & is.null(tsq)){
      tsq <- minque( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs )[1:Q]
  }
  
  if( is.null(theta) ){
    theta <- ginv( t(X)%*%X )%*%(t(X)%*%Y)  
  }
  
  if( null.resids ){
    # Only Anull is needed, so the true constraints are irrelevant
    Anull <- create.constraints( P1, list(order="simple", node=1, decreasing=TRUE) )$Anull
    
    if( is.null(cov.theta) ){
      ssqvec <- rep(ssq,Nks)  
      XSiX   <- t(X) %*% (X/ssqvec)
      if( Q > 0 ){
        tsqvec  <- rep(tsq,Qs)    
        C       <- U * tsqvec
        U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
        tusu       <- t(U1) %*% U1
        diag(tusu) <- diag(tusu) + 1/tsqvec
        tusui      <- solve(tusu)   
        XSiU <- t(X) %*% (U/ssqvec)
        USiY <- t(U) %*% (Y/ssqvec)
        cov.theta <- solve( XSiX - XSiU%*%tusui%*%t(XSiU) )
      } else{
        cov.theta <- solve( XSiX )
      }
    }
    
    if( mySolver=="GLS"){
      wts <- solve(cov.theta)[1:P1, 1:P1, drop=FALSE]
    } else{
      wts <- diag( solve(cov.theta) )[1:P1]
    }
    theta[1:P1] <- activeSet( Anull, y = theta[1:P1], weights = wts, mySolver=mySolver )$x
  }
  
  if( is.null(eps) ){
    resids <- clme_resids( formula, data, gfix )
    eps    <- resids$PA
    xi     <- resids$xi
  }
  
  nu   <- eps
  for( i in 1:K ){
    idx      <- N1[i]:N2[i]
    nu[idx]  <- eps[idx] / sd( eps[idx] )
  }  
  
  ## Obtain the bootstrap samples
  Y.boot  <- matrix( NA, nrow=N, ncol=nsim )
  XT.boot <- X%*%theta
  
  if( Q > 0 ){
    delta <- xi
    for( i in 1:Q ){
      idx         <- Q1[i]:Q2[i]
      delta[idx]  <- xi[idx] / sd( xi[idx] )
    }
    
    Qc <- dim(U)[2]
    for( m in 1:nsim ){
      xi.boot   <- sample( delta, replace=TRUE )
      eps.boot  <- sample( nu   , replace=TRUE )
      
      for( i in 1:Q ){
        idx <- Q1[i]:Q2[i]
        xi.boot[idx]   <- sqrt(tsq[i]) * xi.boot[idx]
      }
      for( i in 1:K ){
        idx <- N1[i]:N2[i]
        eps.boot[idx]  <- sqrt(ssq[i]) * eps.boot[idx]
      }
      Y.boot[,m] <- XT.boot + U%*%xi.boot + eps.boot
    }
    
  } else{
    for( m in 1:nsim ){
      eps.boot  <- sample( nu    , replace=TRUE )
      for( i in 1:K ){
        idx <- N1[i]:N2[i]
        eps.boot[idx]  <- sqrt(ssq[i]) * eps.boot[idx]
      }
      Y.boot[,m] <- XT.boot + eps.boot
    }
  }
  
  ## Return the bootstrap samples
  return( Y.boot )
  
}
