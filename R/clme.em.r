clme_em   <- function( Y, X1, X2 = NULL, U = NULL, Nks = dim(X1)[1],
                     Qs = dim(U)[2], constraints, mq.phi = NULL, tsf = lrt.stat, 
                     tsf.ind = w.stat.ind, mySolver="LS", em.iter = 500, 
                     em.eps =  0.0001, verbose = FALSE, ... ){
  
  
  ##
  ## Development plans: 
  ## - make function inputs more flexible like clme()
  ## - this may make passing the arguments need to be more detailed, but
  ##   passing arguments TO clme.em (e.g. from clme() ) would be less detailed.
  ##
  em_call  <- as.list( environment() )  
  dots     <- as.list(substitute(list(...)))[-1L]
  new_call <- append( em_call, dots )
  
  
  if( is.null(U) ){
    em_results <- do.call( "clme_em_fixed" , new_call )
  } else{
    em_results <- do.call( "clme_em_mixed" , new_call )
  }
  
  return( em_results )

}
