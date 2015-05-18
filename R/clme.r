##
## This is the PACKAGE documentation
##

#' Constrained inference for linear mixed models.
#' 
#' @docType package
#' @name CLME
#' 
#' @description
#' Constrained inference on linear fixed and mixed models using residual bootstrap. 
#' Covariates and random effects are permitted but not required.
#'
#' Appropriate credit should be given when publishing results obtained using \pkg{CLME}, or when 
#' developing other programs/packages based off of this one. Use \code{citation(package="CLME")} for
#'  Bibtex information.
#' 
#' @author Casey M. Jelsema <\email{casey.jelsema@@nih.gov}>
#' @author Shyamal D. Peddada
#' 
#' @note
#' Unless otherwise noted, Casey M. Jelsema wrote the functions in this package.
#' 
#' @details
#' \tabular{ll}{
#' Package: \tab CLME\cr
#' Type:    \tab Package\cr
#' Version: \tab 2.0.3\cr
#' Date:    \tab 2015-05-18\cr
#' License: \tab GLP-3  \cr
#' }
#' 
#' This package implements the constrained linear mixed effects model described in Farnan, Ivanova, and Peddada (2014). See that paper for more details regarding the method. Here we give a brief overview of the assumed model:
#' 
#' \deqn{ Y = X_{1}\theta_{1} + X_{2}\theta_{2} + U\xi + \epsilon }{Y = X1*theta1 + X2*theta2 + U*xi + e}
#' 
#' where
#' 
#' \itemize{
#' \item \eqn{X_1}{X1} is a \eqn{N \times p_1}{N x p1} design matrix.
#' \item \eqn{\theta_1}{theta1} are the coefficients (often treatment effects).
#' \item \eqn{X_2}{X2} is a \eqn{N \times p_2}{N x p2} matrix of fixed covariates.
#' \item \eqn{\theta_1}{theta2} are the coefficients for the covariates.
#' \item \eqn{U}{U} is a \eqn{N \times c}{N x c} matrix of random effects.
#' \item \eqn{\xi}{xi} is a zero-mean random vector with covariance \eqn{T}{T} (see below).
#' \item \eqn{\epsilon}{e} is a zero-mean random vector with covariance \eqn{\Sigma}{Sigma} (see below).
#' }
#' 
#' Neither covariates (\eqn{X_2}{X2}) nor random effects (\eqn{U}{U}) are required by the model or \pkg{CLME}. The covariance matrix of \eqn{\xi}{xi} is given by:
#'   
#' \deqn{ T = diag\left( \tau^{2}_{1}I_{c_{1}}, \tau^{2}_{2}I_{c_{2}}  , \dots , \tau^{2}_{q}I_{c_{q}} \right) }{ T = diag( tau1^2 I_c1, tau2^2 I_c2 , ... , tauq^2 I_cq) }
#' 
#' The first \eqn{c_{1}}{c1} random effects will share a common variance, \eqn{\tau^{2}_{1}}{tau1^2}, the next \eqn{c_{2}}{c2} random effects will share a common variance, and so on. Note that \eqn{c = \sum_{i=1}^{q} c_i}{c = SUM(ci), i=1,...q}. Homogeneity of variances in the random effects can be induced by letting \eqn{q=1}{q=1} (hence \eqn{c_{1}=c=ncol(U)}{c1=c=ncol(U)}).
#' 
#' Similarly, the covariance matrix of \eqn{\epsilon}{e} is given by:
#'   
#' \deqn{ \Sigma = diag\left( \sigma^{2}_{1}I_{n_{1}}, \sigma^{2}_{2}I_{n_{2}}  , \dots , \sigma^{2}_{q}I_{n_{k}} \right) }{ Sigma = diag( sigma1^2 I_n1, sigma2^2 I_n2 , ... , sigmak^2 I_nk)}
#' 
#' Again, the first \eqn{n_{1}}{n1} observations will share a common variance, \eqn{\sigma^{2}_{1}}{sigma1^2}, the next \eqn{n_{2}}{n2} will share a common variance, and so on. Note that \eqn{N = \sum_{i=1}^{k} n_i}{N = SUM(n_i), i=1,...k}. Homogeneity of variances in the residuals can be induced by letting \eqn{k=1}{k=1}.
#' 
#' The order constraints are defined by the matrix \eqn{A}{A}. This is an \eqn{r \times p}{r x p} matrix where \eqn{r}{r} is the number of constraints, and \eqn{p=p_{1}+p_{2}}{p = p1 + p2} is the dimension of \eqn{ \theta = ( \theta_{1}' , \theta_{2}')'}{ theta = ( theta1' , theta2')'}. Formally the hypothesis being tested is:
#'   
#' \deqn{ H_{a}: A\theta > 0 }{Ha: A*theta > 0 }
#' 
#' For several default orders (simple, umbrella, simple tree) the \eqn{A}{A} matrix can be automatically generated. Alternatively, the user may define a custom \eqn{A}{A} matrix to test other patterns among the elements of \eqn{\theta}{theta}. See \code{\link{create.constraints}} and \code{\link{clme}} for more details.
#' 
#' For computational reasons, the implementation is not identical to the model expressed. Particularly, the fixed-effects matrix (or matrices) and the random effects matrix are assumed to be columns in a data frame, not passed as matrices. The \eqn{A}{A} matrix is not \eqn{r\ times p}{r x p}, but \eqn{r\ times 2}{r x 2}, where each row gives the indices of the constrained coefficients. See \code{\link{create.constraints}} for further explanation.
#' 
#' The primary function of \pkg{CLME} is \code{\link{clme}}. The other functions in this package may be run separately, but in general are designed for use by \code{\link{clme}}.
#' 
#' The creation of this package \pkg{CLME}, this manual, and the vingette were all supported by the Intermural Research Program of the United States' National Institutes of Health (Z01 ES101744).
#' 
#' @references
#' Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
#' Linear Mixed Efects Models under Inequality Constraints with Applications.
#' \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
#' \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}
#'
#' @import methods
#' @import stats
#'
#' 
#' 
NULL





#' Constrained Inference for Linear Mixed Effects Models
#'
#' @description Constrained inference for linear fixed or mixed
#'  effects models using distribution-free bootstrap methodology
#'
#' @rdname clme
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
#' @note
#' The argument \code{constraints} is a list containing the order restrictions. The elements are \code{order}, \code{node}, \code{decreasing}, \code{A}, and \code{B}, though not all are necessary. The function can calculate the last two for default orders (simple, umbrella, or simple tree). For default orders, \code{constraints} should be a list containing any subset of \code{order}, \code{node}, and \code{descending}. See the figure below for a depiction of these values; the pictured \code{node} of the simple tree orders (middle column) is 1, and the \code{node} for the umbrella orders (right column) is 3. These may be vectors (e.g. order=('simple','umbrella') ). If any of these three are missing, the function will test for all possible values of the missing element(s), excluding simple tree.
#' For non-default orders, the elements \code{A} and \code{B} should be provided. \code{A} is an \eqn{r \times2}{r x 2} matrix (where r is the number of linear constraints, \eqn{0 < r}{0 < r}. Each row should contain two indices, the first element is the index of the lesser coefficient, the second element is the index of the greater coefficient. So a row of \eqn{(1,2)}{(1,2)} corresponds to the constraint \eqn{\theta_1 \leq \theta_2}{theta_1 <= theta_2}, and a row \eqn{(4,3)}{(4,3)} corresponds to the constraint \eqn{\theta_4 \leq \theta_3}{theta_4 <= theta_3}, etc. Element \code{B} should hold similar contrasts, specifically those needed for calculating the Williams' type test statistic (\code{B} is only needed if \code{tsf=w.stat})
#' The argument \code{tsf} is a function to calculate the desired test statistic. The default function calculates likelihood ratio type test statistic. A Williams type test statistic, which is the maximum of the test statistic over the constraints in \code{constraints\$B}, is also available, and custom functions may be defined. See \code{\link{w.stat}} for details.
#' By default, homogeneity of variances is assumed for residuals (e.g., \code{gfix} does not define groups) and for each random effect.
#' If \code{ncore} is greater than 1 (the default) then \code{clme} will attempt to use the \pkg{doParallel} package to run the bootstrap samples in parallel.
#' 
#' \figure{OrderPlot.jpg}{Plot of Orders.}
#' 
#' @examples
#' data( rat.blood )
#' cons <- list(order="simple", decreasing=FALSE, node=1 )
#' 
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data=rat.blood , 
#'                  constraints=cons, seed=42, nsim=10, ncon=1)
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom MASS ginv
#' 
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
