\name{lrt.stat}

\alias{lrt.stat}


\title{ Likelihood Ratio Type Test Statistic. }
\description{ Calculates the likeihood ratio type test statistic (under Normality assumption) for a constrained linear mixed effects model. This is the default test statistic for \pkg{CLME}. }

\usage{
lrt.stat( theta, theta.null, cov.theta, ... )
}

\arguments{
  \item{theta}{ estimated coefficients.}
  \item{theta.null}{ coefficients estimated under the null hypothesis.}
  \item{cov.theta}{ covariance matrix of the estimated coefficients. For \pkg{CLME} this is taken as the covariance matrix of the unconstrained estimates.}
  \item{\dots}{ additional arguments, to enable custom test statistic functions.}
}


\value{ numeric value }


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This is an internal function, unlikely to be useful outside of \link{CLME-package}. To define custom functions, the arguments available are:

\code{theta}, \code{theta.null}, \code{cov.theta}, \code{B}, \code{A}, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, \code{ssq}, \code{Nks}, and \code{Qs}.

Of the additional arguments, \code{B} and \code{A} are identical to those produced by \code{\link{create.constraints}}. The rest, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, , \code{ssq}, \code{Nks}, and \code{Qs}, are equivalent to arguments to \code{\link{clme_em}}.

Custom functions must produce numeric output. Output may have length greater than 1, which corresponds to testing multiple global hypotheses.

}

\seealso{
\code{\link{clme_em}}, 
\code{\link{w.stat}}
}

\examples{

data( rat.blood )
cons <- list(order = "simple", decreasing = FALSE, node = 1 )

clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
                 constraints = cons, seed = 42, nsim = 0)

# Individually compute lrt statistic
lrt.stat(clme.out$theta, clme.out$theta.null, clme.out$cov.theta )



}



