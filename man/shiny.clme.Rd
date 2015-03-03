\name{shiny_clme}

\alias{shiny_clme}
\alias{shinyServer_clme}
\alias{shinyUI_clme}

\title{ Shiny GUI for CLME }
\description{ A graphical user interface to run \pkg{CLME}, built from the \pkg{shiny} package. }

\usage{
shiny_clme()
shinyServer_clme(input,output)
shinyUI_clme()
}

\arguments{
  \item{input}{ input from GUI.}
  \item{output}{ output to GUI.}
}


\details{
Currently the GUI does not allow specification of custom orders for the alternative hypothesis. Future versions may enable this capability.

The data should be a CSV or table-delimited file with the first row being a header. Variables are identified using their column letter or number (e.g., 1 or A). Separate multiple variables with a comma (e.g., 1,2,4 or A,B,D), or select a range of variables with a dash (e.g., 1-4 or A-D). Set to 'None' (default) to indicate no covariates or random effects.

If group levels for the constrained effect are character, they may not be read in the proper order. An extra column may contain the ordered group levels (it may therefore have different length than the rest of the dataset).

}



\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }



\examples{

\dontrun{ shiny_clme() }

}





