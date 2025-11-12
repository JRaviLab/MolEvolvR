#' InterProScan Column Names
#'
#' A character vector containing the expected column names from an
#' InterProScan output table. This dataset is useful for validating,
#' parsing, or reconstructing data frames produced by InterProScan.
#'
#' @format A character vector with 13 elements:
#' \describe{
#'   \item{AccNum}{Accession number of the sequence.}
#'   \item{SeqMD5Digest}{MD5 digest of the sequence.}
#'   \item{SLength}{Length of the sequence.}
#'   \item{Analysis}{Type of analysis or database used (e.g., Pfam, SMART).}
#'   \item{DB.ID}{Database-specific identifier.}
#'   \item{SignDesc}{Description of the signature or domain.}
#'   \item{StartLoc}{Start position of the match on the sequence.}
#'   \item{StopLoc}{Stop position of the match on the sequence.}
#'   \item{Score}{Score assigned to the match (if applicable).}
#'   \item{Status}{Status of the analysis (e.g., OK, WARNING).}
#'   \item{RunDate}{Date the InterProScan analysis was run.}
#'   \item{IPRAcc}{InterPro accession number.}
#'   \item{IPRDesc}{InterPro entry description.}
#' }
#'
#' @source Generated internally to represent standard InterProScan output fields.
#' @examples
#' data(ipr_colnames)
#' ipr_colnames
"ipr_colnames"
