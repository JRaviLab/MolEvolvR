% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre-msa-tree.R
\name{convert_aln2fa}
\alias{convert_aln2fa}
\title{Adding Leaves to an alignment file w/ accessions}
\usage{
convert_aln2fa(
  aln_file = "",
  lin_file = "data/rawdata_tsv/all_semiclean.txt",
  fa_outpath = "",
  reduced = FALSE
)
}
\arguments{
\item{aln_file}{Character. Path to file. Input tab-delimited file +
alignment file accnum & alignment.
Default is 'pspa_snf7.aln'}

\item{lin_file}{Character. Path to file. Protein file with accession +
number to lineage mapping.
Default is 'pspa.txt'}

\item{fa_outpath}{Character. Path to the written fasta file.
Default is 'NULL'}

\item{reduced}{Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
Default is 'FALSE'}
}
\description{
Adding Leaves to an alignment file w/ accessions
}
\details{
The alignment file would need two columns: 1. accession +
number and 2. alignment. The protein homolog accession to lineage mapping +
file should have
}
\note{
Please refer to the source code if you have alternate +
file formats and/or column names.
}
\examples{
\dontrun{
add_leaves("pspa_snf7.aln", "pspa.txt")
}
}
\author{
Janani Ravi
}
\keyword{accnum,}
\keyword{alignment,}
\keyword{leaves,}
\keyword{lineage,}
\keyword{species}
