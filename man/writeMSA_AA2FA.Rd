% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CHANGED-pre-msa-tree.R, R/pre-msa-tree.R
\name{writeMSA_AA2FA}
\alias{writeMSA_AA2FA}
\title{writeMSA_AA2FA}
\usage{
writeMSA_AA2FA(alignment, outpath)

writeMSA_AA2FA(alignment, outpath)
}
\arguments{
\item{alignment}{MsaAAMultipleAlignment object to be written as a fasta}

\item{outpath}{Where the resulting FASTA file should be written to}
}
\value{
Character string representing the content of the written FASTA file.

Character string of the FASTA content that was written to the file.
}
\description{
MsaAAMultipleAlignment Objects are generated from calls to msaClustalOmega
and msaMuscle from the 'msa' package

Write MsaAAMultpleAlignment Objects as aligned fasta sequence
MsaAAMultipleAlignment Objects are generated from calls to msaClustalOmega
and msaMuscle from the 'msa' package
}
\examples{
\dontrun{
writeMSA_AA2FA("my_sequences.fasta", outpath = "aligned_output.fasta")
}
\dontrun{
# Example usage
alignment <- alignFasta("path/to/sequences.fasta")
writeMSA_AA2FA(alignment, "path/to/aligned_sequences.fasta")
}
}
\author{
Samuel Chen, Janani Ravi
}
