minque <- function( Y , X1 , X2=NULL , U=NULL , Nks=dim(X1)[1] , Qs=dim(U)[2] , 
                    mq.eps=0.0001, mq.iter=500 , verbose=FALSE, ... ){
  
  if( verbose==TRUE ){
    message("Running minque to estimate tau-squared")
  }  
  
  X <- as.matrix( cbind(X1,X2) )
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  K  <- length(Nks)
  
  # Initial values
  theta <- ginv( t(X)%*%X )%*%t(X)%*%Y
  tsq   <- rep( 1 , length(Qs) )
  ssq   <- vector()
  for( k in 1:K ){
    Yk <- Y[ N1[k]:N2[k] ]
    Xk <- X[ N1[k]:N2[k],]  
    ssq[k] <- sum( (Yk-Xk%*%theta)^2 ) / Nks[k]
  }
  
  theta1 <- theta
  tsq1   <- tsq
  ssq1   <- ssq
  
  # Get the F-list
  Flist <- list()
  nullmat <- matrix( 0 , nrow=N , ncol=N )
  if( Q > 0 ){
    for( qk in 1:Q ){
      ind <- Q1[qk]:Q2[qk]
      Flist[[qk]] <- U[,ind] %*% t( U[,ind] )
    }
  }
  
  for( qk in (Q+1):(Q+K) ){
    idx <- N1[qk-Q]:N2[qk-Q]
    X.temp <- nullmat
    diag(X.temp)[idx] <- rep( 1 ,length(idx) )
    Flist[[qk]] <- X.temp  
  }
  
  CONVERGE  <- 0
  iteration <- 0
  phi.hats <- phi.hats1 <- c( tsq , ssq )
  
  # Begin the minque convergence loop
  while( CONVERGE==0 ){
    
    iteration <- iteration + 1
    if( verbose==TRUE ){
      message("Iteration " , iteration )
    }
    
    # Calcualte the G-matrix
    G.phi <- matrix( 0 , nrow=N , ncol=N )
    for( qk in 1:length(Flist) ){
      G.phi <- G.phi + phi.hats1[qk]*Flist[[qk]]
    }
    
    W.phi  <- G.phi + X%*%t(X)
    W.phiI <- ginv(W.phi)
    WX     <- W.phiI%*%X
    R.phi  <- W.phiI - WX%*%ginv( t(X)%*%WX )%*%t(WX)
    YR     <- t(Y)%*%R.phi
    Z.phi  <- as.matrix(sapply( Flist,
                                FUN=function(x, a){ a%*%x%*%t(a) } , a=YR ))
  
    
    S.phi <- matrix( 0 , nrow=(Q+K) , ncol=(Q+K) )
    RFR <- lapply( Flist , FUN=function(x,a){ a%*%x%*%t(a) }, a=R.phi )
    for( qk1 in 1:(Q+K) ){
      S.phi[,qk1] <- sapply( RFR ,
                             FUN=function(x, a){ sum(diag(a%*%x)) } ,
                             a=Flist[[qk1]] )
    }
    
    phi.hats <- ginv(S.phi) %*% Z.phi
    idx.neg  <- phi.hats < 0
    phi.hats[idx.neg,1] <- rep(sqrt(.Machine$double.eps), sum(idx.neg))
    phi.hats <- c(phi.hats)
    
    # Assess convergence    
    rel.change <- abs(phi.hats - phi.hats1)/phi.hats1
    if( mean(rel.change) < mq.eps || iteration >= mq.iter ){
      CONVERGE <- 1
    } else{
      phi.hats1 <- phi.hats  
    }
  
  }
  
  # Return output
  if( Q > 0 ){
    names(phi.hats) <- c(paste("tsq.", 1:Q, sep=""), paste("ssq.", 1:K, sep=""))
  } else{
    names(phi.hats) <- c(paste("ssq.", 1:K,sep=""))    
  }
  
  phi.hats
  
}
