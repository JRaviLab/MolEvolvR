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

#' Map NCBI RefSeq Protein IDs to UniProt Accessions
#'
#' This function maps RefSeq protein accessions (from NCBI) to corresponding
#' UniProt accession IDs using the UniProt ID mapping REST API.
#'
#' @param ncbi_ids A character vector of RefSeq protein accessions (e.g. \code{"NP_001026859.1"}).
#' @param max_retries Integer; maximum number of retries while waiting for the UniProt mapping job (default = 5).
#' @param wait_time Numeric; number of seconds to wait between retries (default = 3).
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{ncbi_id}{Input NCBI RefSeq protein accession(s).}
#'   \item{uniprot_id}{Mapped UniProt accession(s), or \code{NA} if no mapping was found.}
#' }
#'
#' @details
#' This function sends a request to the UniProt ID mapping service
#' (\url{https://rest.uniprot.org/}) and retrieves UniProtKB accessions
#' corresponding to the given RefSeq protein IDs. If UniProt’s service
#' is down or results are delayed, the function retries until \code{max_retries}
#' is reached, after which it returns \code{NA}.
#'
#' @examples
#' \dontrun{
#' ncbi2up(c("NP_001026859.1", "XP_002711597.1"))
#' }
#'
#' @export
ncbi2up <- function(ncbi_ids, max_retries = 5, wait_time = 3) {
    library(httr2)
    library(dplyr)
    library(purrr)
    
    base_url <- "https://rest.uniprot.org/idmapping/run"
    status_base <- "https://rest.uniprot.org/idmapping/status/"
    result_base <- "https://rest.uniprot.org/idmapping/results/"
    
    
    res <- request(base_url) |>
        req_body_form(
            from = "RefSeq_Protein",
            to = "UniProtKB",
            ids = paste(ncbi_ids, collapse = ",")
        ) |>
        req_perform()
    
    res_json <- res |> resp_body_json()
    job_id <- res_json$jobId
    message("Job submitted. Job ID: ", job_id)
    
    
    attempt <- 1
    repeat {
        Sys.sleep(wait_time)
        status <- request(paste0(status_base, job_id)) |>
            req_perform() |>
            resp_body_json()
        
        if (isTRUE(status$jobStatus == "FINISHED") || !is.null(status$results)) break
        
        if (attempt >= max_retries) {
            warning("Timeout: Results not ready. Returning NA.")
            return(tibble(ncbi_id = ncbi_ids, uniprot_id = NA))
        }
        
        message("Waiting for UniProt mapping results... (Attempt ", attempt, ")")
        attempt <- attempt + 1
    }
    
    # Fetching results
    result <- request(paste0(result_base, job_id)) |>
        req_perform() |>
        resp_body_json()
    
    if (is.null(result$results)) {
        warning("No results returned by UniProt API.")
        return(tibble(ncbi_id = ncbi_ids, uniprot_id = NA))
    }
    
    #Extract mappings
    mappings <- map_df(result$results, ~tibble(ncbi_id = .x$from, uniprot_id = .x$to))
    
    
    final <- tibble(ncbi_id = ncbi_ids) |>
        left_join(mappings, by = "ncbi_id")
    
    return(final)
}


#' Map NCBI Protein Accessions to IPG (Identical Protein Group) IDs
#'
#' This function maps NCBI protein accession numbers to their corresponding
#' Identical Protein Group (IPG) IDs using the NCBI Entrez API.
#' It first attempts to retrieve the link via `rentrez::entrez_link()` and,
#' if that fails, falls back to a direct E-utilities REST query.
#'
#' @param acc A character vector of NCBI protein accession numbers.
#' @param api_key (optional) NCBI API key to increase rate limits.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{Accession}{Input NCBI protein accession(s)}
#'   \item{IPG_ID}{Mapped Identical Protein Group ID(s), if available}
#' }
#'
#' @examples
#' \dontrun{
#'   acc2ipg(c("WP_000003915.1", "NP_414543.1"))
#' }
#'
#' @export
acc2ipg <- function(acc, api_key = NULL) {
    
    if (!requireNamespace("rentrez", quietly = TRUE)) {
        stop("Please install 'rentrez' with install.packages('rentrez') or BiocManager::install('rentrez')")
    }
    if (!requireNamespace("httr", quietly = TRUE) ||
        !requireNamespace("xml2", quietly = TRUE)) {
        stop("Please install 'httr' and 'xml2' packages.")
    }
    
    ipg_from_rest <- function(a) {
        url <- paste0(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
            "db=ipg&term=", a, "[Accession]"
        )
        res <- httr::GET(url)
        if (httr::status_code(res) != 200) return(NA)
        xml <- xml2::read_xml(httr::content(res, "text"))
        id <- xml2::xml_text(xml2::xml_find_first(xml, ".//IdList/Id"))
        if (length(id) == 0 || id == "") NA else id
    }
    
    results <- lapply(acc, function(a) {
        tryCatch({
            
            prot_search <- rentrez::entrez_search(
                db = "protein",
                term = paste0(a, "[Accession]"),
                api_key = api_key
            )
            
            if (length(prot_search$ids) == 0) {
                # fallback to REST
                ipg_id <- ipg_from_rest(a)
                return(data.frame(Accession = a, IPG_ID = ipg_id))
            }
            
            
            link_res <- rentrez::entrez_link(
                dbfrom = "protein",
                db = "ipg",
                id = prot_search$ids[1],
                api_key = api_key
            )
            
            ipg_id <- if (!is.null(link_res$links$protein_ipg)) {
                link_res$links$protein_ipg[1]
            } else {
                
                ipg_from_rest(a)
            }
            
            data.frame(Accession = a, IPG_ID = ipg_id)
        },
        error = function(e) {
            message(paste("Error retrieving IPG for", a, ":", e$message))
            data.frame(Accession = a, IPG_ID = NA)
        })
    })
    
    do.call(rbind, results)
}

        
