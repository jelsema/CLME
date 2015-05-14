

clme <-
function( formula, data, gfix=NULL, constraints=list(),
          nsim=1000, tsf=lrt.stat, tsf.ind=w.stat.ind, mySolver="LS", 
          verbose=c(FALSE,FALSE,FALSE), seed=NULL, levels=NULL, ncon=1, ncore=1, ...
          ){
  
  cc       <- match.call( expand.dots=TRUE )  
  
  if( ncon==1 & !is.null(levels) ){
    if( is.list(levels) ){
      idx        <- levels[[1]]
      xlev       <- levels[[2]]
      data[,idx] <- factor( data[,idx] , levels=xlev )
    } else{
      xlev       <- levels
    }
    #idx
  } else{
    xlev <- NULL
  }
  
  mmat     <- model_terms_clme( formula, data, ncon )
  formula2 <- mmat$formula
  Y  <- mmat$Y
  P1 <- mmat$P1
  X1 <- mmat$X1
  X2 <- mmat$X2
  U  <- mmat$U
  
  if( is.null(xlev) ){
    xlev <- colnames(X1)
  } else{
    colnames(X1) <- xlev
  }
  
  if( is.null(gfix) ){
    gfix <- rep("Residual", nrow(X1)) 
  } else{
    data <- with( data, data[order(gfix),])
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
  
  # If only one element for verbose specified, fill the rest with FALSEs
  if( length(verbose)<3 ){
    verbose <- c(verbose, rep(FALSE, 3-length(verbose) ) )
  }  
  
  ## Assess the constraints
  cust.const <- is.matrix( constraints$A )
  prnt_warn <- ""
  
  if( cust.const == TRUE ){
    if( !is.numeric(constraints$A) ){
      stop( "'constraints$A' must be numeric" )
    }    
  } else {
    
    # Constraints are non-null, but A and B are not provided
    # Determine which other elements are missing/needed
    
    if( is.null(constraints$order) ){
      prnt_warn <- paste( prnt_warn, "\n-'constraints$order' is NULL, program will run search for ''simple'' and ''umbrella'' orders")
      constraints$order <- c("simple" , "umbrella" )      
    }
    
    if( is.null(constraints$node) ){
      prnt_warn <- paste( prnt_warn, "\n'constraints$node' is NULL, program will run search for node")
      constraints$node <- 1:P1
    } else{
      search.node <- FALSE
    }
    
    if( is.null(constraints$decreasing) ){
      prnt_warn <- paste( prnt_warn, "\n'constraints$decreasing' is NULL, program will run search for TRUE and FALSE")      
      constraints$decreasing <- c(TRUE,FALSE)      
    }
  }
  
  ## Make sure test stat function is okay
  if( is.function(tsf)==FALSE ){
    stop("'tsf' is not a valid function")
  }
  if( is.function(tsf.ind)==FALSE ){
    stop("'tsf.ind' is not a valid function")
  }
  
  ## Revert to LRT if necessary
  if( cust.const==TRUE & identical( tsf , w.stat ) & is.null(constraints$B) ){
    prnt_warn <- paste( prnt_warn, "\nWilliams type statistic selected with custom constraints, but 
              'constraints$B' is NULL. Reverting to LRT statistic")
    tsf <- lrt.stat
  }
    
  
  ## Set up search grid if using defaults
  if( cust.const==FALSE ){
    search.grid <- expand.grid( constraints$order , 
                                constraints$decreasing ,
                                constraints$node )  
    search.grid[,1] <- as.character(search.grid[,1])
    
    # Remove duplicates / extraneous
    # "simple" doesn't need node
    idx         <- 1*(search.grid[,1]=="simple"   &  search.grid[,3] > 1)
    search.grid <- search.grid[ idx==0 , , drop=FALSE]
    
    # "umbrella" with node=1 or node=P1 is covered by simple order
    if( sum(constraints$order=="simple") >= 1 ){
      idx <- 1*((search.grid[,1]=="umbrella" & search.grid[,3] == 1) + 
                (search.grid[,1]=="umbrella" & search.grid[,3] == P1))
      search.grid <- search.grid[ idx==0 , , drop=FALSE]
    } else{
      idx <- 1*(search.grid[,1]=="umbrella" & search.grid[,3] == 1)
      search.grid[idx,1] <- rep( "simple" , sum(idx) )
      idx <- 1*(search.grid[,1]=="umbrella" & search.grid[,3] == P1)
      search.grid <- search.grid[ idx==0 , , drop=FALSE]
    }
    
    # Move simple.tree to the bottom
    idx <- search.grid[,1]=="simple.tree"
    if( sum(idx)>0 ){
      search.grid <- rbind( search.grid[idx==0, , drop=FALSE] ,
                            search.grid[idx==1, , drop=FALSE] )
    }
    
    ## A check for duplicate rows here may be wise
    MNK <- dim( search.grid )[1]  
    
  } else{
    MNK <- 1
    loop.const <- est.const <- constraints
  }
  
  
  ##
  ## End preparation steps, begin the analysis
  ##
  
  ## Obtain tau if needed
  if( is.null(U)==FALSE ){
    mq.phi <- minque( Y=Y , X1=X1 , X2=X2 , U=U , Nks=Nks , Qs=Qs ,
                      verbose=verbose[2], ... )
  } else{
    mq.phi <- NULL
  }
  
  ## EM for the observed data
  if( verbose[1]==TRUE ){
    print( paste( "Starting EM Algorithm for observed data." , sep=""))
  }
  
  ## Loop through the search grid
  est.order <- NULL
  ts.max    <- -Inf

  
  for( mnk in 1:MNK ){
    
    if( cust.const==FALSE ){
      grid.row <- list( order     = search.grid[mnk,1], 
                        node      = search.grid[mnk,3], 
                        decreasing= search.grid[mnk,2])       
      loop.const <- create.constraints( P1=ncol(X1), constraints=grid.row  )
    }
        
    clme.temp <- clme_em( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, 
                          Qs=Qs, constraints=loop.const, mq.phi=mq.phi,
                          tsf=tsf, tsf.ind=tsf.ind, mySolver=mySolver,
                          verbose=verbose[3], ... )    
    
    # If global test stat is larger, update current estimate of order  
    if( cust.const==FALSE ){
      update.max <- (mnk==1) + (clme.temp$ts.glb > ts.max)
    } else{
      update.max <- 1
    }
    
    if( update.max > 0 ){
      ts.max    <- clme.temp$ts.glb
      clme.out  <- clme.temp
      est.order <- mnk
    }
    
  }
  
  if( cust.const==FALSE ){
    grid.row <- list( order     = search.grid[est.order,1], 
                      node      = search.grid[est.order,3],
                      decreasing= search.grid[est.order,2]) 
    est.const <- create.constraints( P1=ncol(X1), constraints=grid.row  ) 
  } else{
    est.const <- constraints
  }
  
  ## Calculate the residuals from unconstrained model
  mr <- clme_resids( formula=formula, data=data, gfix=gfix, ncon=ncon )
  
  ## This is the loop for the bootstrap simulations
  if( nsim > 0 ){
    if( round(ncore)==ncore & ncore > 1   ){
      
      registerDoParallel( cores=ncore )
      
      ## Use PARALLEL processing for the bootstrap simulations
      ## Use SEQUENTIAL processing for the bootstrap simulations
      ## Obtain bootstrap samples      
      Y.boot <- resid_boot( formula=formula, data=data, gfix=gfix, 
                            eps=mr$PA, xi=mr$xi, ssq=mr$ssq, tsq=mr$tsq, 
                            cov.theta=mr$cov.theta, nsim=nsim, 
                            theta=clme.out$theta.null, mySolver=mySolver,
                            seed=seed, null.resids=FALSE, ncon=ncon, ...  )
      
      ## EM for the bootstrap samples    
      #p.value  <- rep( 0 , length(clme.out$ts.glb) )
      #pval.ind <- rep( 0 , dim(est.const$A)[1] )
      p.value  <- matrix( 0, nrow=nsim , ncol= )
      pval.ind <- matrix( 0, nrow=nsim , ncol=dim(est.const$A)[1]     )
      
      mprint <- round( seq( 1 , round(nsim*0.9), length.out=10 ) )
      
      
      pvals <- foreach( m = 1:nsim , .combine='rbind' ) %dopar% {
        
        if( verbose[1]==TRUE & (m %in% mprint) ){
          print( paste( "Bootstrap Iteration " , m , " of " , nsim , sep=""))
        }
        
        ## Loop through the search grid
        ts.boot <- -Inf
        
        for( mnk in 1:MNK ){
          if( cust.const==FALSE ){
            grid.row <- list( order=search.grid[mnk,1], node=search.grid[mnk,3],
                              decreasing=search.grid[mnk,2] )
            loop.const <- create.constraints( P1=ncol(X1), constraints=grid.row )
          }
          
          clme.temp <- clme_em( Y=Y.boot[,m], X1=X1, X2=X2, U=U, Nks=Nks,
                                Qs=Qs, constraints=loop.const, mq.phi=mq.phi,
                                tsf=tsf, tsf.ind=tsf.ind, mySolver=mySolver,
                                verbose=verbose[3], ...)
          
          idx <- which(clme.temp$ts.glb > ts.boot)
          if( length(idx)>0 ){
            ts.boot[idx] <- clme.temp$ts.glb[idx]
          }
          
          update.ind <- (MNK==1) + (mnk == est.order)
          if( update.ind>0 ){
            ts.ind.boot <- clme.temp$ts.ind 
          }
        }
        c( 1*( ts.boot >= clme.out$ts.glb ), 1*(ts.ind.boot >= clme.out$ts.ind) )
      }
      
      p.value  <- pvals[ ,   1:length(clme.out$ts.glb) , drop=FALSE ]
      pval.ind <- pvals[ , -(1:length(clme.out$ts.glb)), drop=FALSE ]
      
      clme.out$p.value     <- colSums(p.value)/nsim
      clme.out$p.value.ind <- colSums(pval.ind)/nsim
      
      ## End of the PARALLEL BOOTSTRAP LOOP
    } else{
      ## Use SEQUENTIAL processing for the bootstrap simulations
      ## Obtain bootstrap samples      
      Y.boot <- resid_boot( formula=formula, data=data, gfix=gfix, 
                            eps=mr$PA, xi=mr$xi, ssq=mr$ssq, tsq=mr$tsq, 
                            cov.theta=mr$cov.theta, nsim=nsim, 
                            theta=clme.out$theta.null, mySolver=mySolver,
                            seed=seed, null.resids=FALSE, ncon=ncon, ...  )
      
      ## EM for the bootstrap samples    
      p.value  <- rep( 0 , length(clme.out$ts.glb) )
      pval.ind <- rep( 0 , dim(est.const$A)[1] )
      
      mprint <- round( seq( 1 , round(nsim*0.9), length.out=10 ) )
      
      for( m in 1:nsim ){
        
        if( verbose[1]==TRUE & (m %in% mprint) ){
          print( paste( "Bootstrap Iteration " , m , " of " , nsim , sep=""))
        }
        
        ## Loop through the search grid
        ts.boot <- -Inf
        
        for( mnk in 1:MNK ){
          if( cust.const==FALSE ){
            grid.row <- list( order=search.grid[mnk,1], node=search.grid[mnk,3],
                              decreasing=search.grid[mnk,2] )
            loop.const <- create.constraints( P1=ncol(X1), constraints=grid.row )
          }
          
          clme.temp <- clme_em( Y=Y.boot[,m], X1=X1, X2=X2, U=U, Nks=Nks,
                                Qs=Qs, constraints=loop.const, mq.phi=mq.phi,
                                tsf=tsf, tsf.ind=tsf.ind, mySolver=mySolver,
                                verbose=verbose[3], ...)
          
          idx <- which(clme.temp$ts.glb > ts.boot)
          if( length(idx)>0 ){
            ts.boot[idx] <- clme.temp$ts.glb[idx]
          }
          
          update.ind <- (MNK==1) + (mnk == est.order)
          if( update.ind>0 ){
            ts.ind.boot <- clme.temp$ts.ind 
          }
        }
        p.value  <- p.value  + 1*( ts.boot    >= clme.out$ts.glb )
        pval.ind <- pval.ind + 1*(ts.ind.boot >= clme.out$ts.ind )
      }
      
      clme.out$p.value     <- p.value/nsim
      clme.out$p.value.ind <- pval.ind/nsim
      
    } ## End of the SEQUENTIAL BOOTSTRAP LOOP
  } else{
    clme.out$p.value     <- NA
    clme.out$p.value.ind <- rep( NA, nrow(est.const$A) )
  }

  
  ## Add some values to the output object
  class(clme.out)       <- "clme"
  clme.out$constraints  <- list( A=est.const$A, B=est.const$B )
  clme.out$dframe        <- mmat$dframe
  
  names(clme.out$theta) <- c( colnames(X1), colnames(X2) )
  names(clme.out$ssq)   <- names(Nks)
  names(clme.out$tsq)   <- names(Qs)
  
  if( !is.null(levels) ){
    names(clme.out$theta)[1:P1]        <- xlev
    colnames(clme.out$cov.theta)[1:P1] <- xlev
    rownames(clme.out$cov.theta)[1:P1] <- xlev
  }
  
  if( is.null(U) ){
    clme.out$residuals    <- mr$PA
  } else{
    clme.out$residuals    <- cbind( mr$PA, mr$SS, mr$FM )
    colnames(clme.out$residuals) <- c("PA", "SS", "FM")
  }
  
  clme.out$random.effects <- mr$xi

  clme.out$gfix    <- Nks
  clme.out$gran    <- Qs
  clme.out$formula <- mmat$formula
  clme.out$call    <- cc  
  clme.out$P1      <- P1
  clme.out$nsim    <- nsim
  
  ## Report the estimated order
  clme.out$order <- list()
  if( cust.const == TRUE ){
    clme.out$order$estimated <- FALSE
    clme.out$order$order     <- "custom"
    clme.out$order$node      <- NULL
    clme.out$order$inc.dec   <- NULL
  } else{
      if( MNK==1 ){
        clme.out$order$estimated <- FALSE
      } else{
        clme.out$order$estimated <- TRUE
      }
      
      clme.out$order$order <- est.const$order
      clme.out$order$node  <- est.const$node  
      
      if( est.const$decreasing ){
        clme.out$order$inc.dec <- "decreasing"
      } else{
        clme.out$order$inc.dec <- "increasing"
      }
      
      
  }
  
  if (verbose[1]==TRUE){
    cat( prnt_warn )
  }
  
  ## Return the output object
  return( clme.out )
  
}
