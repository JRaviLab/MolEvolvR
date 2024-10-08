% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.R
\name{words2WordCounts}
\alias{words2WordCounts}
\title{Words 2 Word Counts}
\usage{
words2WordCounts(string)
}
\arguments{
\item{string}{}
}
\value{
\link{tbl_df} table with 2 columns: 1) words & 2) counts/frequency
}
\description{
Get word counts (wc) \link{DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)}
}
\examples{
\dontrun{
tibble::tibble(DomArch = c("aaa+bbb", "a+b", "b+c", "b-c")) |>
    elements2Words() |>
    words2WordCounts()
}

}
