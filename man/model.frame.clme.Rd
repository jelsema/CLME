\name{model.frame.clme}

\alias{model.frame.clme}


\title{ Extracts the model frame }
\description{
Extracts the model frame from objects of class \code{clme}. 
}

\usage{
\method{model.frame}{clme}(formula, ... )
}

\arguments{
  \item{formula}{ a formula expression. }
  \item{...}{space for additional arguments (should include argument \code{data}).}  
}


\value{
Returns a data frame with the variables in the model.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{clme}}
}

\examples{

\dontrun{
data( rat.blood )
model.frame.clme( mcv ~ time + temp + sex + (1|id), data = rat.blood )
}

}

