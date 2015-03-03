

##
## Williams' type statistic (global)
##
w.stat <- function( theta , cov.theta , B , A , ...  ){
  
  stats <- vector( "numeric",length=nrow(B) )
  ctd   <- diag( cov.theta )
  
  stats <- apply( B , 1 , 
                  FUN=function(b,theta,cov,ctd){
                    std <- sqrt( ctd[b[1]] + ctd[b[2]] - 2*cov.theta[b[1],b[2]] )
                    (theta[b[2]]-theta[b[1]])/std
                  }, theta=theta, cov=cov.theta, ctd=ctd)
  
  test.stat <- max( stats )
  tnames    <- names(theta)
  ii        <- min( which( stats == max(stats) ) )  
  
  if( !is.null(tnames) ){
    names(test.stat) <- paste0( tnames[B[ii,2]] , " - ", tnames[B[ii,1]] )
  } else{
    names(test.stat) <- paste0( "Theta ", B[ii,2] , " - Theta ", B[ii,1] )
  }
  
  return( test.stat )
  
}


##
## Williams' type statistic (individual)
##
w.stat.ind <- function( theta , cov.theta , B , A , ...  ){
  
  stats <- vector( "numeric",length=nrow(A) )
  ctd   <- diag( cov.theta )
  
  stats <- apply( A , 1 , 
                  FUN=function(a,theta,cov,ctd){
                    std <- sqrt( ctd[a[1]] + ctd[a[2]] - 2*cov.theta[a[1],a[2]] )
                    (theta[a[2]]-theta[a[1]])/std
                  }, theta=theta, cov=cov.theta, ctd=ctd)
  
  return( stats )
  
}


##
## Likelihood ratio type statistic (global)
##
lrt.stat <- function( theta, theta.null, cov.theta, ... ){  
  theta.diff <- theta - theta.null  
  test.stat  <- c( t(theta.diff) %*% cov.theta %*% theta.diff )
  
  names(test.stat) <- "Bootstrap LRT"
  
  # Return test statistic
  return(test.stat)
  
}


