#' These functions will be deprecated. Please use other functions instead.
#' 
#' @name deprecate
#' 
NULL

#' @rdname deprecate
#' @export
sink.reset <- function() {
    warning("'sink.reset' is deprecated. Use 'sinkReset' instead.")
    sinkReset() 
}

#' @rdname deprecate
#' @export
add_lins <- function(df, ...) {
    warning("'add_lins' is deprecated. Use 'addlineage' instead.")
    addlineage(df, ...) 
}

#' @rdname deprecate
#' @export
acc2lin <- function(accessions, ...) {
    warning("'acc2lin' is deprecated. Use 'acc2Lineage' instead.")
    acc2Lineage(accessions, ...)
}

#' @rdname deprecate
#' @export
efetch_ipg <- function(accnums, ...) {
    warning("'efetch_ipg' is deprecated. Use 'efetchIPG' instead.")
    efetchIPG(accnums, ...) 
}

#' @rdname deprecate
#' @export
ipg2lin <- function(accessions, ...) {
    warning("'ipg2lin' is deprecated. Use 'IPG32Lineage' instead.")
    IPG32Lineage(accessions, ...)  
}