create.constraints <- function( P1, constraints ){
  
  Q1 <- P1-1

  
  order      <- tolower(constraints$order)
  node       <- constraints$node
  decreasing <- constraints$decreasing
  
  if( is.null(node) ){
    node <- 1
  }
  
  if( (order %in% c("simple", "simple.tree", "umbrella"))==FALSE ){
    stop("'order' must be one or more of: simple, simple.tree, umbrella")
  }
  
  ## Revert to simple order if umbrella has node at extreme
  if( order=="umbrella" & node %in% c(1,P1) ){
    order <- "simple"
    node  <- 1
    if( node==P1 ){
      if( decreasing==TRUE ){
        decreasing <- FALSE
      } else {
        decreasing <- TRUE
      }
    }
  }
  
  
  A <- matrix( 0, nrow=Q1, ncol=2 )
  
  ## Simple order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_K
  if( order=="simple" ){
    if( decreasing==TRUE ){
      A <- as.matrix(cbind( 1:Q1+1 , 1:Q1 ))
      B <- matrix( c(P1,1) , nrow=1 )
    } else{
      A <- as.matrix(cbind( 1:Q1 , 1:Q1 + 1 ))
      B <- matrix( c(1,P1) , nrow=1 )
    }
    node    <-  NULL
  }
  
  ## Simple tree order
  ## e.g. mu_1 <= mu_i ; i=2,...,P1
  if( order=="simple.tree" ){
    if( decreasing==TRUE ){
      A <- as.matrix( cbind( (1:P1)[-node] , rep(node,Q1) ) )
    } else{
      A <- as.matrix(cbind( rep(node,Q1)  , (1:P1)[-node] ))
    }
    B <- A    
  }
  
  ## Umbrella order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_b >= mu_{b+1} >= ... >= mu_K
  if( order=="umbrella" ){    
    if( decreasing==TRUE ){
      for( ii in 1:(node-1) ){
        A[ii,] <- c(ii,ii+1)
      }
      for( ii in node:Q1 ){
        A[ii,] <- c(ii+1,ii)
      }
      B <- as.matrix( rbind( c(1,node) , c(P1,node) ) )
      
    } else{
      
      for( ii in 1:(node-1) ){
        A[ii,] <- c(ii+1,ii)
      }
      for( ii in node:Q1 ){
        A[ii,] <- c(ii,ii+1)
      }
      B <- as.matrix( rbind( c(node,1) , c(node,P1) ) )      
    }
  }
  
  
  ## Make the null A-matrix
  Anull <- matrix( 0 , nrow=P1*(P1-1) , ncol=2 )
  for( ii in 1:P1 ){
    idx1 <- (P1-1)*(ii-1) + 1
    idx2 <- (P1-1)*ii
    Anull[idx1:idx2,1] <- rep(ii,P1-1)
    Anull[idx1:idx2,2] <- (1:P1)[-ii]
  }
  
    
  
  # Return the constraints object
  new_constraints <- list( A = A, B = B, Anull = Anull, order=order, node=node, decreasing=decreasing )
  return(new_constraints)
}
