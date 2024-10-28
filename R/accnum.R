

#' Convert UniProt IDs to NCBI RefSeq Accessions
#'
#' This function takes one or more UniProt IDs and returns the corresponding NCBI RefSeq accessions.
#' It uses the org.Hs.eg.db package to perform the mapping.
#'
#' @param uniprot_ids A character vector of UniProt IDs.
#' @return A data frame with columns 'UNIPROT' and 'REFSEQ', mapping UniProt IDs to RefSeq accessions.
#'         Returns an empty data frame if no mappings are found.
#' @examples
#' \dontrun{
#'   uniprot_ids <- c("P04217", "P01023")
#'   refseq_accessions <- up2ncbi(uniprot_ids)
#'   print(refseq_accessions)
#' }
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @export
up2ncbi <- function(uniprot_ids) {
  # Check if input is provided
  if (length(uniprot_ids) == 0) {
    stop("No UniProt IDs provided.")
  }
  
  # Perform the mapping
  tryCatch({
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = uniprot_ids,
      columns = "REFSEQ",
      keytype = "UNIPROT"
    )
    
    # Check if any mappings were found
    if (nrow(mapping) == 0) {
      warning("No NCBI accessions found for the given UniProt IDs.")
    }
    
    return(mapping)
  }, error = function(e) {
    stop(paste("Error in mapping:", e$message))
  })
}



#' Convert NCBI RefSeq Accessions to UniProt IDs
#'
#' This function takes one or more NCBI RefSeq accession numbers and returns the corresponding UniProt IDs.
#' It uses the org.Hs.eg.db package to perform the mapping.
#'
#' @param ncbi_accessions A character vector of NCBI RefSeq accession numbers.
#' @return A data frame with columns 'REFSEQ' and 'UNIPROT', mapping RefSeq accessions to UniProt IDs.
#'         Returns NA for UniProt if no mapping is found.
#' @examples
#' \dontrun{
#'   ncbi_accessions <- c("NP_000005.2", "NP_000007.1")
#'   uniprot_ids <- ncbi2up(ncbi_accessions)
#'   print(uniprot_ids)
#' }
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @export
ncbi2up <- function(ncbi_accessions) {
  # Check if input is provided
  if (length(ncbi_accessions) == 0) {
    stop("No NCBI accessions provided.")
  }
  
  # Strip version numbers from accessions
  stripped_accessions <- gsub("\\.[0-9]+$", "", ncbi_accessions)
  
  # Perform the mapping
  tryCatch({
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = stripped_accessions,
      columns = "UNIPROT",
      keytype = "REFSEQ"
    )
    
    # Check if any mappings were found
    if (nrow(mapping) == 0) {
      warning("No UniProt IDs found for the given NCBI accessions.")
    }
    
    return(mapping)
  }, error = function(e) {
    stop(paste("Error in mapping:", e$message))
  })
}
