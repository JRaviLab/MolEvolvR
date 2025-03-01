% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotStackedLineage}
\alias{plotStackedLineage}
\title{Stacked Lineage Plot}
\usage{
plotStackedLineage(
  prot,
  column = "DomArch",
  cutoff,
  Lineage_col = "Lineage",
  xlabel = "Domain Architecture",
  reduce_lineage = TRUE,
  label.size = 8,
  legend.position = c(0.7, 0.4),
  legend.text.size = 10,
  legend.cols = 2,
  legend.size = 0.7,
  coord_flip = TRUE,
  legend = TRUE,
  cpcols = NULL
)
}
\arguments{
\item{prot}{Data frame containing protein data including domain architecture
and lineage information.}

\item{column}{Character. The name of the column in prot representing domain
architectures (default is "DomArch").}

\item{cutoff}{Numeric. A threshold value for filtering domain architectures
or protein counts.}

\item{Lineage_col}{Character. The name of the column representing lineage
data (default is "Lineage").}

\item{xlabel}{Character. Label for the x-axis
(default is "Domain Architecture").}

\item{reduce_lineage}{Logical. Whether to shorten lineage names
(default is TRUE).}

\item{label.size}{Numeric. The size of axis text labels (default is 8).}

\item{legend.position}{Numeric vector. Coordinates for placing the legend
(default is c(0.7, 0.4)).}

\item{legend.text.size}{Numeric. Size of the text in the legend
(default is 10).}

\item{legend.cols}{Numeric. Number of columns in the legend (default is 2).}

\item{legend.size}{Numeric. Size of the legend keys (default is 0.7).}

\item{coord_flip}{Logical. Whether to flip the coordinates of the plot
(default is TRUE).}

\item{legend}{Logical. Whether to display the legend (default is TRUE).}

\item{cpcols}{}
}
\value{
A ggplot object representing a stacked bar plot showing the
distribution of protein domain architectures across lineages.
}
\description{
Stacked Lineage Plot
}
\examples{
\dontrun{
plotStackedLineage()
}
}
