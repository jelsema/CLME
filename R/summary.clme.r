#' S3 method to summarize results for objects of class \code{clme}
#'
#' @description Summarizes the output of objects of class \code{clme}, such as those produced by \code{\link{clme}}. Prints a tabulated display of global and individual tests, as well as parameter estimates.
#'
#' @param object an object of class \code{clme}.
#' @param alpha level of significance.
#' @param digits number of decimal digits to print.
#' @param ... additional arguments passed to other functions.
#'
#' @note 
#' The individual tests are performed on the specified order. If no specific order was specified, then the individual tests are performed on the estimated order.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#'   set.seed( 42 )
#'   data( rat.blood )
#'   cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#'   clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                    constraints = cons, seed = 42, nsim = 10)
#'   
#'   summary( clme.out )
#' }
#' 
#' @importFrom stringr str_pad
#' @importFrom stringr str_trim
#' @importFrom prettyR decimal.align
#' 
#' @method summary clme
#' @export
#' 
summary.clme <- function( object, alpha=0.05, digits=4, ...){
  
  ## Title and formula
  cat( "Linear mixed model subject to order restrictions\n" )
  cat( "Formula: ")
  print( object$formula )
  
  if( object$order$order=="unconstrained" ){
    cat( paste0("\nNo order restrictions (two-tailed alternatives)") )
  } else{
    ## Order statement
    if( object$order$order=="simple" ){
      order <- "simple order"
    }
    if( object$order$order=="umbrella" ){
      order <- paste0("umbrella order with node at ", object$order$node )
    }
    if( object$order$order=="simple.tree" ){
      order <- paste0("tree order with node at ", object$order$node )
    }
    
    if( object$order$order == "custom" ){
      cat( "\nCustom order constraints were provided" )
    } else{
      if( object$order$estimated ){
        ## Estimated the order
        cat( paste0("\nOrder estimated: " , object$order$inc.dec , " ", order ) )
      } else{
        cat( paste0("\nOrder specified: " , object$order$inc.dec , " ", order ) )
      }
    }
  }
  
    
  ## Diagnostic criterion
  crit <- c(logLik.clme(object),
            AIC.clme(object),
            AIC.clme( object, k=log(nobs.clme(object)/(2*pi)) ) )
  critc <- format( crit , digits=4)
  cat( "\n\nlog-likelihood:", critc[1] )
  cat( "\nAIC:           "  , critc[2] )
  cat( "\nBIC:           "  , critc[3] )
  cat( "\n(log-likelihood, AIC, BIC computed under normality)")  
  
  ## Tests
  est    <- fixef.clme(object)
  tnames <- names(est)
  Amat   <- object$constraints$A
  Bmat   <- object$constraints$B

  ## Global tests
  if( object$order$order != "unconstrained" ){
  if( length(object$ts.glb)>1 ){
    glbs <- object$ts.glb
    grow <- matrix( "NA" , nrow=length(glbs), ncol=3 )
    
    for( ii in 1:length(glbs) ){
      #if( is.null(Bmat) ){
      #  glbn <- "Unknown"
      #  glbe <- "Unknown"
      #} else{
      #  glbn <- paste( tnames[Bmat[ii,2]] , "-", tnames[Bmat[ii,1]] , sep=" " )
      #  glbe <- round( est[Bmat[ii,2]] - est[Bmat[ii,1]], digits=3 )
      #}
      ## Estimate will generally be NA anyway (for LRT test), just drop it.
      # grow[ii,] <- c( glbn, glbe, round(object$ts.glb[ii],3) , sprintf("%.4f", object$p.value[ii]) )
      if( is.null(names(glbs)) ){
        glbn <- rep( "Unknown" , length(glbs) )
      } else{
        glbn <- names( object$ts.glb )
      }
      grow[ii,] <- c( glbn[ii], round(object$ts.glb[ii],3) , sprintf("%.4f", object$p.value[ii]) )
      
    }
    
    #colnames( grow ) <- c("Contrast", "Estimate", "Stat", "p-value")
    #grow <- .align_table.clme( grow )
    for( ii in 2:3){
      val1 <- str_trim( decimal.align( grow[,ii]), side="right"  )
      grow[,ii] <- str_pad(val1, width=max(nchar(val1)), side = "right", pad = "0")
    }  
    
    #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
    colnames( grow ) <- c("Contrast", "Statistic", "p-value")
    grow1 <- c(colnames(grow)[1], grow[,1])
    grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
    grow2 <- .align_table.clme( grow[,2:3,drop=FALSE] )
    grow <- cbind( grow1[2:length(grow1)] , grow2)
    colnames(grow)[1] <- grow1[1]
    
    cat( "\n\nGlobal tests: ")
    cat( "\n", paste(colnames(grow) , collapse="  ") )
    for( ii in 1:length(glbs) ){
      cat( "\n", paste(grow[ii,] , collapse="  ")     )
    }
    
  } else{
    ## Single global tests
    #if( is.null(Bmat) ){
    #  glbn <- "Unknown"
    #  glbe <- "Unknown"
    #} else{
    #  glbn <- paste( tnames[Bmat[1,2]] , "-", tnames[Bmat[1,1]] )
    #  glbe <- round( est[Bmat[1,2]] - est[Bmat[1,1]], digits=3 )
    #}
    
    #grow <- cbind( glbn, glbe, round(object$ts.glb,3) , sprintf("%.4f", object$p.value) )
    #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
    if( is.null(names(object$ts.glb)) ){
      glbn <- "Unknown"
    } else{
      glbn <- names( object$ts.glb )
    }
    grow <- cbind( glbn, round(object$ts.glb,3) , sprintf("%.4f", object$p.value) )
    
    colnames( grow ) <- c("Contrast", "Statistic", "p-value")
    grow1 <- c(colnames(grow)[1], grow[,1])
    grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
    grow2 <- .align_table.clme( grow[,2:3,drop=FALSE] )
    grow <- cbind( grow1[2:length(grow1)] , grow2)
    colnames(grow)[1] <- grow1[1]
    
    cat( "\n\nGlobal test: ")
    cat( "\n", paste0(colnames(grow) , collapse="  ") )   
    cat( "\n", paste0(grow , collapse="  ")     )
  }
  }
  
  ## Individual tests
  glbs <- object$ts.ind
  grow <- matrix( "NA" , nrow=length(glbs), ncol=4 )
  
  for( ii in 1:length(glbs) ){
    glbn <- paste( tnames[Amat[ii,2]] , "-", tnames[Amat[ii,1]] )
    glbe <- round( est[Amat[ii,2]] - est[Amat[ii,1]], digits=3 )
    grow[ii,] <- c( glbn, glbe, round(object$ts.ind[ii],3) , sprintf("%.4f", object$p.value.ind[ii]) )
  }
  
  #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
  #grow <- .align_table.clme( grow[,2:4] )

  for( ii in 2:4){
    
    if( !any(grow[,ii]==rep("NA",nrow(grow))) ){
      val1 <- str_trim( decimal.align( grow[,ii]), side="right"  )
      grow[,ii] <- str_pad(val1, width=max(nchar(val1)), side = "right", pad = "0")   
    }
  }  
  
  colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
  grow1 <- c(colnames(grow)[1], grow[,1])
  grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
  grow2 <- .align_table.clme( grow[,2:4,drop=FALSE] )
  grow <- cbind( grow1[2:length(grow1)] , grow2)
  colnames(grow)[1] <- grow1[1]
    
  cat( "\n\nIndividual Tests (Williams' type tests): ")
  cat( "\n", paste(colnames(grow) , collapse="  ") )
  for( ii in 1:length(glbs) ){
    cat( "\n", paste(grow[ii,] , collapse="  ")     )
  }

  
  ## Random effects
  cat( "\n\nVariance components: \n")
  print( VarCorr.clme(object) )
  
  ## Fixed effects
  vars  <- diag( vcov.clme(object) )
  CIs   <- confint.clme( object , level=(1-alpha), ...)
  tvals <- cbind(tnames = str_pad(tnames, width=max(nchar(tnames)), side = "right", pad = " "), 
                 cest   = format(est       , digits=4),
                 cvars  = format(sqrt(vars), digits=4),
                 clcl   = format(CIs[,1]   , digits=4),
                 cucl   = format(CIs[,2]   , digits=4))
    
  
  cipct <- round(100*(1-alpha),2)
  colnames(tvals) <- c(" ", "Estimate", "Std. Err",
                       paste0(cipct, "% lower"), paste0(cipct, "% upper"))
  tvals <- .align_table.clme( tvals )
  
  cat( "\n\nFixed effect coefficients (theta): \n")
  cat( paste0(colnames(tvals), collapse="  ") )
  for( ii in 1:length(est) ){
    cat( "\n", paste0( c(tvals[ii,]),  collapse="  ") )
  }
  cat( "\nStd. Errors and confidence limits based on unconstrained covariance matrix")
  
  cat( "\n\nParameters are ordered according to the following factor levels:\n" )
  cat( paste(  names(fixef(object))[1:object$P1], collapse=", ") )
  cat( "\n\nModel based on", paste0(object$nsim), "bootstrap samples" )
  
}


