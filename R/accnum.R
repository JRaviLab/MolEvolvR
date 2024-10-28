

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



#' Convert NCBI Protein Accessions to IPG (Identical Protein Group) IDs
#'
#' This function takes one or more NCBI protein accession numbers and returns the corresponding
#' IPG (Identical Protein Group) IDs. It uses the NCBI E-utilities API via the rentrez package
#' to perform the mapping.
#'
#' @param ncbi_ids A character vector of NCBI protein accession numbers.
#' @return A data frame with columns 'NCBI' and 'IPG', mapping NCBI protein accessions to IPG IDs.
#'         Returns an empty data frame if no mappings are found.
#' @examples
#' \dontrun{
#'   ncbi_ids <- c("NP_000005.2", "NP_000007.1")
#'   ipg_mappings <- ncbi2ipg(ncbi_ids)
#'   print(ipg_mappings)
#' }
#' @importFrom rentrez entrez_search
#' @export
ncbi2ipg <- function(ncbi_ids) {
  # Check if input is provided
  if (length(ncbi_ids) == 0) {
    stop("No NCBI IDs provided.")
  }
  
  # Perform the mapping
  tryCatch({
    # Search the IPG database for each NCBI ID
    results <- lapply(ncbi_ids, function(id) {
      search <- entrez_search(db = "ipg", term = paste0(id, "[PACC]"))
      if (search$count > 0) {
        data.frame(NCBI = id, IPG = search$ids, stringsAsFactors = FALSE)
      } else {
        NULL
      }
    })
    
    # Combine results into a single data frame
    mapping <- do.call(rbind, results)
    
    # Check if any mappings were found
    if (is.null(mapping) || nrow(mapping) == 0) {
      warning("No IPG mappings found for the given NCBI IDs.")
      mapping <- data.frame(NCBI = character(0), IPG = character(0), stringsAsFactors = FALSE)
    }
    
    return(mapping)
  }, error = function(e) {
    stop(paste("Error in mapping:", e$message))
  })
}