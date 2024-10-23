



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