#' Map UniProt IDs to NCBI accessions
#'
#' This function queries UniProt and returns the corresponding NCBI protein accessions.
#' If the UniProt API is unreachable, it returns NA values with a warning.
#'
#' @param uniprot_ids A character vector of UniProt accession IDs.
#' @return A tibble with UniProt IDs and mapped NCBI accessions.
#' @examples
#' up2ncbi_fast(c("P12345", "Q9Y263"))
#' @export
up2ncbi <- function(uniprot_ids) {
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:(", paste(uniprot_ids, collapse = " OR "), ")")

    resp <- tryCatch(
        {
            httr2::request(base_url) |>
                httr2::req_url_query(query = query, format = "json", fields = "accession,xref_refseq") |>
                httr2::req_perform()
        },
        error = function(e) {
            warning("UniProt API request failed: ", conditionMessage(e))
            return(NULL)
        }
    )

    if (is.null(resp)) {
        return(tibble::tibble(uniprot_id = uniprot_ids,
                              ncbi_accession = NA_character_))
    }

    results <- httr2::resp_body_json(resp)

    if (length(results$results) == 0) {
        return(tibble::tibble(uniprot_id = uniprot_ids,
                              ncbi_accession = NA_character_))
    }

    df <- tibble::tibble(
        uniprot_id = vapply(results$results, function(x) x$primaryAccession, ""),
        ncbi_accession = vapply(results$results, function(x) {
            if (!is.null(x$uniProtKBCrossReferences)) {
                paste0(
                    unique(unlist(lapply(x$uniProtKBCrossReferences, function(y) {
                        if (y$database == "RefSeq") y$id else NULL
                    }))),
                    collapse = ";"
                )
            } else NA_character_
        }, "")
    )

    df <- dplyr::right_join(df, tibble::tibble(uniprot_id = uniprot_ids), by = "uniprot_id")

    return(df)
}
