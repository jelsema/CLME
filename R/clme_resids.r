clme_resids <-
function( formula, data, gfix=NULL, ncon=1 ){
  
  ##
  ## I should consider a way to condense this so it doesn't need to 
  ## be copied between here, clme() and resid_boot().
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
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix( cbind(X1, X2) )
  K  <- length(Nks)
  
  
  # Initial values
  theta <- ginv( t(X)%*%X )%*%( t(X)%*%Y )
  
    ssq   <- vector()
    for( k in 1:K ){
      Yk <- Y[ N1[k]:N2[k] ]
      Xk <- X[ N1[k]:N2[k],]  
      ssq[k] <- sum( (Yk - Xk%*%theta)^2 ) / (Nks[k])
    }
  
  ## Obtain the estimates of epsilon and delta
  ssqvec  <- rep(ssq,Nks)    
  XSiX <- t(X) %*% (X/ssqvec)
  XSiY <- t(X) %*% (Y/ssqvec)
  
  if( Q > 0 ){
    mq.phi <- minque( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs )[1:Q]
    tsq    <- mq.phi
    
    tsqvec  <- rep(tsq,Qs)    
    C       <- U * tsqvec
    
    U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve(tusu)   
    XSiU <- t(X) %*% (U/ssqvec)
    USiY <- t(U) %*% (Y/ssqvec)
    XPiX <- XSiX - XSiU%*%tusui%*%t(XSiU)
    XPiY <- XSiY - XSiU%*%(tusui%*%USiY)
  } else{
    XPiX <- XSiX
    XPiY <- XSiY
  }
  
  # H <- X%*%ginv( XPiX )%*%(t(X)%*%PsiI)
  #eps  <- c( Y - H%*%Y )  
  Yhat <- X %*% ginv( XPiX )%*%XPiY
  PA   <- c(Y - Yhat)
  
  if( Q > 0 ){
    # xi    <- c( t(C) %*% PsiI %*% eps )  
    USiR <- t(U) %*% (PA/ssqvec)
    USiU <- t(U) %*% (U/ssqvec)
    xi   <- tsqvec * ( USiR - USiU%*%(tusui%*%USiR) )
    SS <- U %*% xi
    FM <- PA - SS

    resid.out <- list( PA=c(PA), SS=c(SS), FM=c(FM), cov.theta=solve(XPiX),
                       xi=c(xi), ssq=ssq, tsq=tsq )
  } else{
    resid.out <- list( PA=c(PA), cov.theta=XPiX, ssq=ssq )
  }
  
  return( resid.out )
    
}












