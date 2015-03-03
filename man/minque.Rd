\name{minque}
\alias{minque}
\title{ MINQUE Algorithm }
\description{ Algorithm to obtain MINQUE estimates of variance components of a linear mixed effects model.}

\usage{
minque( Y , X1 , X2 = NULL , U = NULL , Nks = dim(X1)[1] ,
        Qs = dim(U)[2], mq.eps =  0.0001 ,
        mq.iter = 500 , verbose = FALSE, ... )
}

\arguments{
  \item{Y}{ \eqn{N \times 1}{N x 1}  vector of response data. }
  \item{X1}{ \eqn{N \times p_1}{N x p1} design matrix. }
  \item{X2}{ optional \eqn{N \times p_2}{N x p2} matrix of covariates. }
  \item{U}{ optional \eqn{N \times c}{N x c}  matrix of random effects. }
  \item{Nks}{ optional \eqn{K \times 1}{K x 1} vector of group sizes. See Details. }
  \item{Qs}{ optional \eqn{Q \times 1}{Q x 1} vector of group sizes for random effects.}
  \item{mq.eps}{ criterion for convergence for the MINQUE algorithm. }
  \item{mq.iter}{ maximum number of iterations permitted for the MINQUE algorithm. }
  \item{verbose}{ if \code{TRUE}, function prints intermediate messages on progress of the MINQUE algorithm. }
  \item{...}{space for additional arguments.}
}

\details{
By default, the model assumes homogeneity of variances for both the residuals and the random effects (if included). See the Details in \code{\link{clme_em}} for more information on how to use the arguments \code{Nks} and \code{Qs} to permit heterogeneous variances. 
}

\value{
The function returns a vector of the form \eqn{(\tau^{2}_{1}, \tau^{2}_{2}, \ldots, \tau^{2}_{q}, \sigma^{2}_{1},\sigma^{2}_{2},\ldots, \sigma^{2}_{k})'}{(tau1^2, tau2^2, \ldots, tauq^2, sigma1^2,sigma2^2,\ldots, sigmak^2)'}. If there are no random effects, then the output is just \eqn{(\sigma^{2}_{1},\sigma^{2}_{2},\ldots, \sigma^{2}_{k})'}{(sigma1^2,sigma2^2,\ldots, sigmak^2)'}.
}


\references{
Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}
}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This function is called by several other function in \pkg{CLME} to obtain estimates of the random effect variances. If there are no random effects, they will not call \code{minque}.
}

\examples{

data( rat.blood )

model_mats <- model_terms_clme( mcv ~ time + temp + sex + (1|id) , 
                                data = rat.blood )

Y  <- model_mats$Y
X1 <- model_mats$X1
X2 <- model_mats$X2
U  <- model_mats$U

# No covariates or random effects
minque(Y = Y, X1 = X1 )

# Include covariates and random effects
minque(Y = Y, X1 = X1, X2 = X2, U = U )


}


