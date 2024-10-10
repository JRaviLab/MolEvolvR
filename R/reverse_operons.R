# Function to straighten operons (genomic contexts)
# Written by L. Aravind
# Modified by Janani Ravi and Samuel Chen


#' reveql: Reverse Equalities in Genomic Context
#'
#' @description
#' This function processes the genomic context strings (GenContext) and reverses
#'  directional signs based on the presence of an equal sign ("="). 
#'
#' @param prot [vector] A vector of genomic context strings to be processed.
#'
#' @return [vector] A vector of the same length as the input, where each genomic 
#' element is annotated with either a forward ("->") or reverse ("<-") direction, 
#' depending on its position relative to the "=" symbols.
#'
#' @export
#'
#' @examples
#' # Example input: Genomic context with directional symbols and an asterisk
#' genomic_context <- c("A", "B", "*", "C", "D", "=", "E", "F")
#' reveql(genomic_context)
#'
#' # Output: "A->", "B->", "*", "<-C", "<-D", "=", "E->", "F->"
reveql <- function(prot) {
    w <- prot # $GenContext.orig # was 'x'

    y <- rep(NA, length(w))

    d <- 1

    b <- grep("\\*", w)

    for (j in b:length(w)) {
        if (w[j] == "=") {
            d <- d * (-1)
        }

        if (d == 1 && w[j] != "=") {
            y[j] <- paste(w[j], "->", sep = "")
        } else if (d == -1 && w[j] != "=") {
            y[j] <- paste("<-", w[j], sep = "")
        } else {
            y[j] <- "="
        }
    } # (for)

    if (b > 1) {
        d <- 1

        for (j in (b - 1):1) {
            if (w[j] == "=") {
                d <- d * (-1)
            }

            if (d == 1 && w[j] != "=") {
                y[j] <- paste(w[j], "->", sep = "")
            } else if (d == -1 && w[j] != "=") {
                y[j] <- paste("<-", w[j], sep = "")
            } else {
                y[j] <- "="
            }
        } # (for)
    } # (if b>1)

    return(y)
}

## The function to reverse operons

#' reverse_operon: Reverse the Direction of Operons in Genomic Context
#'
#' @description
#' This function processes a genomic context data frame to reverse the direction
#' of operons based on specific patterns in the GenContext column. It handles 
#' elements represented by ">" and "<" and restructures the genomic context by 
#' flipping the direction of operons while preserving the relationships 
#' indicated by "=".
#'
#' @param prot [data.frame] A data frame containing at least a column named 
#' 'GenContext', which represents the genomic contexts that need to be reversed.
#'
#' @return [data.frame] The input data frame with the 'GenContext' column updated t
#' o reflect the reversed operons.
#'
#' @export
#'
#' @examples
#' # Example genomic context data frame
#' prot <- data.frame(GenContext = c("A>B", "C<D", "E=F*G", "H>I"))
#' reversed_prot <- reverse_operon(prot)
#' print(reversed_prot)
reverse_operon <- function(prot) {
    gencontext <- prot$GenContext

    gencontext <- gsub(pattern = ">", replacement = ">|", x = gencontext)

    gencontext <- gsub(pattern = "<", replacement = "|<", x = gencontext)

    gencontext <- gsub(pattern = "\\|\\|", replacement = "\\|=\\|", x = gencontext)



    gc.list <- strsplit(x = gencontext, split = "\\|")

    if (any(is.na(gc.list))) gc.list[[which(is.na(gc.list))]] <- "-"

    gc.list <- lapply(1:length(gc.list), function(x) {
        if (any(gc.list[[x]] == "")) gc.list[[x]][which(gc.list[[x]] != "")] else gc.list[[x]]
    })



    te <- lapply(1:length(gc.list), function(x) gc.list[[x]][grep("\\*", gc.list[[x]])])

    ye <- unlist(lapply(te, function(x) substr(x[1], 1, 1)))

    torev <- which(ye == "<")



    te <- gc.list[torev]

    te <- lapply(te, function(x) gsub(pattern = "<-|->", replacement = "", x = x))

    te <- lapply(te, rev)

    witheq <- grep(pattern = "=", x = te)

    withouteq <- which(!((1:length(te)) %in% witheq))

    ge <- te[witheq]



    ge <- lapply(1:length(ge), function(x) reveql(ge[[x]]))

    ye <- te[withouteq]

    ye <- lapply(1:length(ye), function(x) unname(sapply(ye[[x]], function(y) paste(y, "->", sep = ""))))



    te[witheq] <- ge

    te[withouteq] <- ye

    gc.list[torev] <- te



    rev.gencontext <- unlist(lapply(gc.list, function(x) paste(x, collapse = "")))

    rev.gencontext <- gsub(pattern = "=", replacement = "\\|\\|", rev.gencontext)

    prot$GenContext <- rev.gencontext

    return(prot)
}



##############
# Absorb into the function above?
## ???
# colnames(prot) <- c("AccNum","GenContext.orig","len", "GeneName","TaxID","Species")

## ??? straighten operons
# prot$GenContext.orig <- reverse_operon(prot)
