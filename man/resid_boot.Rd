\name{resid_boot}

\alias{resid_boot}

\title{Obtain Residual Bootstrap}

\description{ Generates bootstrap samples of the data vector. }


\usage{resid_boot(formula, data, gfix=NULL, eps=NULL, xi=NULL, null.resids=TRUE,
                  theta=NULL, ssq=NULL, tsq=NULL, cov.theta=NULL, seed=NULL, 
                  nsim=1000, mySolver="LS", ncon=1, ... ) }

\arguments{
  \item{formula}{ a formula expression. The constrained effect(s) must come before any unconstrained covariates on the right-hand side of the expression. The first \code{ncon} terms will be assumed to be constrained.}
  \item{data}{ data frame containing the variables in the model. }
  \item{gfix}{ optional vector of group levels for residual variances. Data should be sorted by this value.}
  \item{eps}{ estimates of residuals. }
  \item{xi}{ estimates of random effects. }
  \item{null.resids}{ logical indicating if residuals should be computed under the null hypothesis. }
  \item{theta}{estimates of fixed effects coefficients. Estimated if not submitted.}
  \item{ssq}{estimates of residual variance components. Estimated if not submitted.}
  \item{tsq}{estimates of random effects variance components. Estimated if not submitted.}
  \item{cov.theta}{covariance matrix of fixed effects coefficients. Estimated if not submitted.}
  \item{seed}{ set the seed for the RNG. }
  \item{nsim}{ number of bootstrap samples to generate. }
  \item{mySolver}{solver to use, passed to \code{activeSet}.}
  \item{ncon}{ the number of variables in \code{formula} that are constrained.}
  \item{...}{space for additional arguments.}
}


\details{
If any of the parameters \code{theta}, \code{ssq}, \code{tsq}, \code{eps}, or \code{xi} are provided, the function will use those values in generating the bootstrap samples. They will be estimated if not submitted. If\code{null.resids=TRUE}, then \code{theta} will be projected onto the space of the null hypothesis ( \eqn{H_{0}: \theta_1 = \theta_2 = ... = \theta_{p_1}}{Ho: theta_1 = theta_2 = ... = theta_p1}) regardless of whether it is provided or estimated. To generate bootstraps with a specific \code{theta}, set \code{null.residuals=FALSE}.


}


\value{ Output is \eqn{N \ times nsim}{N x nsim} matrix, where each column is a bootstrap sample of the response data \code{Y}. }

\references{
Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}
}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This function is primarily designed to be called by \code{\link{clme}}. 

By default, homogeneous variances are assumed for the residuals and (if included) random effects. Heterogeneity can be induced using the arguments \code{Nks} and \code{Qs}, which refer to the vectors \eqn{ (n_{1}, n_{2}, \ldots, n_{k}) }{(n1, n2 ,... , nk)} and \eqn{ (c_{1}, c_{2}, \ldots, c_{q}) }{(c1, c2 ,... , cq)}, respectively. See \code{\link{clme_em}} for further explanation of these values.

}

\seealso{
\code{\link{clme}}
}


\examples{

data( rat.blood )
boot_sample <- resid_boot(mcv ~ time + temp + sex + (1|id), nsim = 10, 
                         data = rat.blood, null.resids = TRUE, ncon = 1 )
                        
                        
                        

}


