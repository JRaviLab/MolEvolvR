# Function to straighten operons (genomic contexts)
# Written by L. Aravind
# Modified by Janani Ravi and Samuel Chen


#' straightenOperonSeq
#'
#' @param prot
#'
#' @return
#' @export
#'
#' @examples
straightenOperonSeq <- function(prot) {
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

#' reverseOperonSeq
#'
#' @param prot
#'
#' @return
#' @export
#'
#' @examples
reverseOperonSeq <- function(prot) {
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



    ge <- lapply(1:length(ge), function(x) straightenOperonSeq(ge[[x]]))

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
# prot$GenContext.orig <- reverseOperonSeq(prot)
