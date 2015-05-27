
#' Unconstrained Inference for Linear Mixed Effects Models
#'
#' @description Linear fixed or mixed
#'  effects models using distribution-free bootstrap methodology
#'
#' @rdname pw_clme
#'
#' @param formula a formula expression. The constrained effect(s) must come before any unconstrained covariates on the right-hand side of the expression. The first \code{ncon} terms will be assumed to be constrained.
#' @param data data frame containing the variables in the model. 
#' @param gfix optional vector of group levels for residual variances. Data should be sorted by this value.
#' @param constraints optional list containing the constraints. See Details for further information. 
#' @param nsim optional number of bootstrap samples to use for significance testing. 
#' @param tsf function to calculate the test statistic. 
#' @param tsf.ind function to calculate the test statistic for individual constrats. See Details for further information. 
#' @param mySolver solver to use in isotonization (passed to \code{activeSet}). 
#' @param verbose optional. Vector of 3 logicals. The first causes printing of iteration step, the second two are passed as the \code{verbose} argument to the functions \code{\link{minque}} and \code{\link{clme_em}}, respectively. 
#' @param seed set the seed for the RNG.
#' @param levels optional list to manually specify names for constrained coefficients. See Details.
#' @param ncon the number of variables in \code{formula} that are constrained.
#' @param ncore the number of cores to use in parallel processing.
#' @param ... space for additional arguments.
#'
#'
#' @details 
#' If any random effects are included, the function computes MINQUE estimates of variance components. After, \code{\link{clme_em}} is run to obtain the observed values. If \code{nsim}>0, a bootstrap test is performed using \code{\link{resid_boot}}.
#' For the argument \code{levels} the first list element should be the column index (in \code{data}) of the constrained effect. The second element should be the true order of the levels.
#'
#' @return
#' The output of \code{clme} is an object of the class \code{clme}, which is list with elements:
#' \itemize{    
#' \item{\code{theta}}{ estimates of \eqn{\theta}{theta} coefficients}
#' \item{\code{theta}}{ estimates of \eqn{\theta_0}{theta_0} coefficients under the null hypothesis}
#' \item{\code{ssq}}{ estimate of residual variance(s), \eqn{\sigma^{2}_{i}}{sigma.i^2}.}
#' \item{\code{tsq}}{ estimate of random effects variance component(s), \eqn{\tau^{2}_{i}}{tau.i^2}.}
#' \item{\code{cov.theta}}{ the unconstrained covariance matrix of \eqn{\theta}{theta}}
#' \item{\code{ts.glb}}{ test statistic for the global hypothesis.}
#' \item{\code{ts.ind}}{ test statistics for each of the constraints.}
#' \item{\code{mySolver}}{ the solver used for isotonization.}
#' \item{\code{p.value}}{ p-value for the global hypothesis}
#' \item{\code{p.value.ind}}{ p-values for each of the constraints}
#' \item{\code{constraints}}{ list containing the constraints (\code{A}) and the contrast for the global test (\code{B}).}
#' \item{\code{dframe}}{ data frame containing the variables in the model.}
#' \item{\code{residuals}}{ matrix containing residuals. For mixed models three types of residuals are given. }
#' \item{\code{random.effects}}{ estimates of random effects. }
#' \item{\code{gfix}}{ group sample sizes for residual variances. }
#' \item{\code{gran}}{ group sizes for random effect variance components. }
#' \item{\code{formula}}{ the formula used in the model. }
#' \item{\code{call}}{ the function call. }
#' \item{\code{order}}{ list describing the specified constraints.}
#' \item{\code{P1}}{ the number of constrained parameters.}
#' \item{\code{nsim}}{ the number of bootstrap simulations used for inference.}
#' }
#' 
#' 
#' @examples
#' data( rat.blood )
#' cons <- list(order="simple", decreasing=FALSE, node=1 )
#' 
#' clme.out <- pw_clme(mcv ~ time + temp + sex + (1|id), data=rat.blood , 
#'                  constraints=cons, seed=42, nsim=10, ncon=1)
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom MASS ginv
#' 
#' @export
#' 
pw_clme <-
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
  
  ## Make sure test stat function is okay
  if( is.function(tsf)==FALSE ){
    stop("'tsf' is not a valid function")
  }
  if( is.function(tsf.ind)==FALSE ){
    stop("'tsf.ind' is not a valid function")
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
  
  pairwise <- list( A     = t( combn(1:ncol(X1), m=2) ) , 
                    Anull = create.constraints( P1, constraints=list(order="simple",node=1,decreasing=TRUE) )$Anull)
  
  clme.out <- pw_clme_em( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, 
                        Qs=Qs, constraints=pairwise, mq.phi=mq.phi, mySolver=mySolver,
                        tsf=tsf, tsf.ind=tsf.ind, verbose=verbose[3], ... )    
  
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
      mprint <- round( seq( 1 , round(nsim*0.9), length.out=10 ) )
      
      pvals <- foreach( m = 1:nsim , .combine='rbind', .packages="foreach" ) %dopar% {
        
        if( verbose[1]==TRUE & (m %in% mprint) ){
          print( paste( "Bootstrap Iteration " , m , " of " , nsim , sep=""))
        }
        
        ## Loop through the search grid        
        clme.temp <- pw_clme_em( Y=Y.boot[,m], X1=X1, X2=X2, U=U, Nks=Nks,
                              Qs=Qs, constraints=pairwise, mq.phi=mq.phi,
                              tsf=tsf, tsf.ind=tsf.ind, mySolver=mySolver,
                              verbose=verbose[3], ...)
        
        ts.glb <- clme.temp$ts.glb
        ts.ind <- clme.temp$ts.ind
        
        pvglb <- 1*( (ts.glb >= abs(clme.out$ts.glb) ) | (ts.glb <= -abs(clme.out$ts.glb) )  )
        pvind <- 1*( (ts.ind >= abs(clme.out$ts.ind) ) | (ts.ind <= -abs(clme.out$ts.ind) )  )
                
        c( pvglb, pvind )
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
      pval.ind <- rep( 0 , nrow(pairwise$A) )
      
      mprint <- round( seq( 1 , round(nsim*0.9), length.out=10 ) )
      
      for( m in 1:nsim ){
        
        if( verbose[1]==TRUE & (m %in% mprint) ){
          print( paste( "Bootstrap Iteration " , m , " of " , nsim , sep=""))
        }
                
        clme.temp <- pw_clme_em( Y=Y.boot[,m], X1=X1, X2=X2, U=U, Nks=Nks,
                              Qs=Qs, constraints=pairwise, mq.phi=mq.phi,
                              tsf=tsf, tsf.ind=tsf.ind, mySolver=mySolver,
                              verbose=verbose[3], ...)
        
        ts.glb <- clme.temp$ts.glb
        ts.ind <- clme.temp$ts.ind
        
        p.value  <- p.value + 1*( (ts.glb >= abs(clme.out$ts.glb) ) | (ts.glb <= -abs(clme.out$ts.glb) )  )
        pval.ind <- pval.ind + 1*( (ts.ind >= abs(clme.out$ts.ind) ) | (ts.ind <= -abs(clme.out$ts.ind) )  )
        
      }
      
      clme.out$p.value     <- p.value/nsim
      clme.out$p.value.ind <- pval.ind/nsim
      
    } ## End of the SEQUENTIAL BOOTSTRAP LOOP
  } else{
    clme.out$p.value     <- NA
    clme.out$p.value.ind <- rep( NA, nrow(pairwise$A) )
  }

  
  ## Add some values to the output object
  class(clme.out)       <- "clme"
  clme.out$constraints  <- list( A=pairwise$A, B=NULL )
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
  clme.out$order$order     <- "unconstrained"
  
  ## Return the output object
  return( clme.out )
  
}
