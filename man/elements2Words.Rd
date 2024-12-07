% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.R
\name{elements2Words}
\alias{elements2Words}
\title{elements2Words}
\usage{
elements2Words(prot, column = "DomArch", conversion_type = "da2doms")
}
\arguments{
\item{prot}{A dataframe containing the dataset to analyze. The specified
\code{column} contains the string elements to be processed.}

\item{column}{A character string specifying the name of the column to analyze.
Default is "DomArch".}

\item{conversion_type}{A character string specifying the type of conversion.
Two options are available:
\describe{
\item{\code{da2doms}}{Convert domain architectures into individual domains by
replacing \code{+} symbols with spaces.}
\item{\code{gc2da}}{Convert genomic context into domain architectures by
replacing directional symbols (\verb{<-}, \verb{->}, and \code{|}) with spaces.}
}}
}
\value{
A single string where elements are delimited by spaces. The function
performs necessary substitutions based on the \code{conversion_type} and cleans up
extraneous characters like newlines, tabs, and multiple spaces.
}
\description{
Break string ELEMENTS into WORDS for domain architecture (DA) and genomic
context (GC)
}
\examples{
\dontrun{
tibble::tibble(DomArch = c("aaa+bbb", 
"a+b", "b+c", "b-c")) |> elements2Words()
}

}
