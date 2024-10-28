



#' Convert a UniProt ID to an NCBI Entrez Gene ID
#'
#' This function takes a single UniProt ID and returns the corresponding NCBI Entrez Gene ID.
#' It uses the `org.Hs.eg.db` package to perform the mapping.
#' 
#' @author Klangina
#' @param uniprot_id A string representing a single UniProt ID.
#' @return A string representing the corresponding NCBI Entrez Gene ID. Returns `NA` if no mapping is found.
#' @examples
#' \dontrun{
#'   uniprot_id <- "P04217"
#'   entrez_id <- up2ncbi(uniprot_id)
#'   print(entrez_id)
#' }
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @export
up2ncbi <- function(uniprot_id) {
  # Use the select function to map the UniProt ID to an Entrez Gene ID
  result <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = uniprot_id,
                                  columns = "ENTREZID",
                                  keytype = "UNIPROT")
  
  # Check if the result is not empty and return the first Entrez ID
  if (nrow(result) > 0 && !is.na(result$ENTREZID[1])) {
    return(as.character(result$ENTREZID[1]))
  } else {
    return(NA)  # Return NA if no mapping is found
  }
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
