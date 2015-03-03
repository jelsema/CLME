



plot.clme <- 
  function(x , alpha=0.05 , legendx="below" , inset=0.01,
           ci=FALSE , ylim=NULL , cex=1.75 , pch=21 , bg="white" , 
           xlab = expression( paste( "Component of " , theta[1] ) ),
           ylab = expression( paste( "Estimated Value of " , theta[1] ) ) , 
           tree=NULL, ...){
  object <- x
  if( !is.clme(object) ){
    stop("Argument 'object' is not of class clme.")
  } else{
    
    theta <- fixef(object)
    A     <- object$constraints$A
    r     <- nrow(A)
    p1    <- max(A)
    
    if( ci ){ 
      ci.wd <- min( 1/p1 , 1/15 ) 
      CIs <- confint(object, level=(1-alpha))
    }
    
    if( legendx=="below" ){
      layout( rbind(1,2) , heights=c(7,1) )
    }
    
    theta1 <- object$theta[1:p1]
    
    # Pick some reasonable plot limits
    if( is.null(ylim) ){
      if( ci ){
        ylim <- c( min(CIs[1:p1,]) , max(CIs[1:p1,]))
      } else{
        if( min(theta1) < 0 ){ 
          ymin <- min(theta1)*1.05 
        } else{
          ymin <- min(theta1)/1.05
        }
        if( max(theta1) < 0 ){ 
          ymax <- max(theta1)/1.05
        } else{
          ymax <- max(theta1)*1.05
        }
        ylim <- c( ymin , ymax )
      }
    }
    
    if( is.null(tree) ){
      tree <- object$order$order == "simple.tree"      
    }
    
    
    ## PLOT FOR SIMPLE / UMBRELLA ORDERS
    if( !tree ){
      # The initial plot of the points
      plot( 1:p1 , theta1 , cex=cex , pch=pch , bg=bg ,
          ylim = ylim , xaxt='n' , xlab = xlab , ylab = ylab, ...)
      
      axis(side=1, at=1:p1, labels=names(theta1), ...)
  
      # Connect the contrasts with solid/dashed lines
      for( ii in 1:r){  
        idx <- A[ii,]
        if( object$p.value.ind[ii]  > alpha ){ lty <- 1 }
        if( object$p.value.ind[ii] <= alpha ){ lty <- 2 }
        points( idx , theta[idx] , lty=lty , lwd=2 , type="l")
      }
  
      ## Add the CIs if necessary
      if( ci ){
        for( ii in 1:p1){
          points( c(ii,ii)  , c(CIs[ii,1], CIs[ii,2]) , type="l"  )
          points( c(ii-ci.wd, ii+ci.wd)  , c(CIs[ii,1], CIs[ii,1]) , type="l"  )
          points( c(ii-ci.wd, ii+ci.wd)  , c(CIs[ii,2], CIs[ii,2]) , type="l"  )        
        }    
      }
      
      # Replot the pointsso the circles are filled
      points( 1:p1 , theta[1:p1] , cex=cex , pch=pch , bg=bg )
    }
    
    if( tree ){
      ## MAKE THE CODE FOR A SIMPLE TREE PLOT HERE.      
      plot(x=1, y=0, col=0, ylim=ylim, xlim=c(0.9,2.1), xlab="", ylab="Estimated Coefficient", xaxt="n")
      axis(side=1, at=c(1,1.78), labels=c("Control (Node)" , "Treatment") )
      node <- object$order$node
      
      legend( 0.86, theta1[node]+0.35, names(theta1)[node] ,cex=.8, bty='n' ) 
      
      for( ii in (1:p1)[-node] ){
          legend(   1.77    , theta1[ii]+0.15, names(theta1)[ii] ,cex=.8, bty='n' ) 
          points( c(1,1.78) , theta1[c(node,ii)] , col=1 , type="l" , lwd=2 , lty=(1 + 1*(object$p.value.ind[ii-1] < alpha)) )
          points( c(1,1.78) , theta1[c(node,ii)] , col=1 , cex=1.5 , pch=21 , bg="white" )        
      }
    }
    
    
    ## Put a legend on the plot if requested
    if( legendx=="below" ){
      safe.mar <- par( no.readonly=TRUE )$mar
      par(mar=c(0, 0, 0, 0) , ...)
      plot.new()
      legend('center','groups', c( paste("p >" , alpha, "  " ) , paste("p <" , alpha, "  " )),
             lty = c(1,2), col=1 , ncol=2 , bty ="o" , ...)
      
      par( mar=safe.mar )
    } else{
      leg.texts <- c("bottom", "bottomleft", "left", "topleft",
                    "top", "topright", "right", "bottomright", "center")
      if( legendx %in% leg.texts){
        legend( legendx , legend=c( paste("p >" , alpha, "  " ) , paste("p <" , alpha, "  " )), 
                lty = c(1,2), col=1 , ncol=1 , bty ="o" , inset=inset , ...)
      }
    }
    
  }
  
}


