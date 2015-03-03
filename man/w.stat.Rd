\name{w.stat}

\alias{w.stat}
\alias{w.stat.ind}


\title{ Williams' Type Test Statistic. }
\description{ Calculates a Williams' type test statistic for a constrained linear mixed effects model. }

\usage{
w.stat( theta , cov.theta , B , A , ...)
w.stat.ind( theta , cov.theta , B , A , ...)
}

\arguments{
  \item{theta}{ estimated coefficients.}
  \item{cov.theta}{ covariance matrix of the (unconstrained) coefficients. }
  \item{B}{   matrix to obtain the global contrast. }
  \item{A}{   matrix of linear constraints.}
  \item{\dots}{ additional arguments, to enable custom test statistic functions.}
}


\details{
See \code{\link{create.constraints}} for an example of \code{A}. Argument \code{B} is similar, but defines the global contrast for a Williams' type test statistic. This is the largest hypothesized difference in the constrained coefficients. So for an increasing simple order, the test statistic is the difference between the two extreme coefficients, \eqn{\theta_1}{theta_1} and \eqn{\theta_{p_1}}{theta_p1}, divided by the standard error (unconstrained). For an umbrella order order, two contrasts are considered, \eqn{\theta_1}{theta_1} to \eqn{\theta_{s}}{theta_s}, and \eqn{\theta_{p_1}}{theta_p1} to \eqn{\theta_{s}}{theta_s}, each divided by the appropriate unconstrained standard error. A general way to express this statistic is:

\deqn{W = max \theta_{B[i,2]} - \theta_{B[i,1]} / sqrt( VAR( \theta_{B[i,2]} - \theta_{B[i,1]} ) )}{W = max theta_{B[i,2]} - theta_{B[i,1]} / sqrt( VAR( theta_{B[i,2]} - theta_{B[i,1]} ) )}

where the numerator is the difference in the constrained estimates, and the standard error in the denominator is based on the covariance matrix of the unconstrained estimates.

The function \code{w.stat.ind} does the same, but uses the \code{A} matrix which defines all of the individual constraints, and returns a test statistic for each constraints instead of taking the maximum.


%If \eqn{B}{B} is a matrix  with each row \eqn{B_i}{Bi} being a contrast, and \eqn{\Sigma}{Sig} is the covariance matrix of \eqn{\theta}{theta}, \code{w.stat} returns:
%\deqn{W = max \frac{B_i \theta}{ B_{i} \Sigma^{-1} B_{i}^{'} }}{ W = max [Bi * theta] / [ Bi * Sigma^-1 * Bi' }
%and \code{w.stat.ind} returns the vector W where:
%\deqn{W_{i} = \frac{A_i \theta}{ A_{i} \Sigma^{-1} A_{i}^{'} }}{ Wi = max [Ai * theta] / [ Ai * Sigma^-1 * Ai' }

}



\value{ Output is a numeric value. }
\references{

Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}

}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
See \code{\link{lrt.stat}} for information on creating custom test statistics.
}


\examples{

theta  <- exp(1:4/4)
th.cov <- diag(4)
X1     <- matrix( 0 , nrow=1 , ncol=4 )
const  <- create.constraints( P1=4 , constraints=list(order='simple' ,
          decreasing=FALSE) )

w.stat( theta , th.cov , const$B , const$A )

w.stat.ind( theta , th.cov , const$B , const$A )

}






