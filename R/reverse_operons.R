#' Reverse Operon Sequences in Genomic Contexts
#'
#' This function takes a data frame with genomic contexts and processes them
#' to reverse the direction of operons in sequences containing specific
#' symbols such as `>`, `<`, and `=`.
#'
#' @param prot A data frame with a column named 'GenContext' containing genomic contexts.
#' @return A modified data frame with the operons' directions reversed.
#' @export
#'
#' @examples
#' # Example genomic context data frame
#' prot <- data.frame(GenContext = c("A>B", "C<D", "E=F*G", "H>I"))
#' reversed_prot <- reverseOperonSeq(prot)
#' print(reversed_prot)
reverseOperonSeq <- function(prot) {
  # Ensure input is in correct format (GenContext must be a character vector)
  if (!"GenContext" %in% colnames(prot) || !is.character(prot$GenContext)) {
    stop("GenContext must be a character column in the input data frame")
  }

  # Modify GenContext to ensure proper splitting and replacements
  gencontext <- prot$GenContext
  gencontext <- gsub(">", ">|", gencontext)
  gencontext <- gsub("<", "|<", gencontext)
  gencontext <- gsub("\\|\\|", "|=|", gencontext)

  # Split GenContext into a list and remove empty elements
  gc.list <- strsplit(gencontext, "\\|")
  gc.list <- lapply(gc.list, function(x) x[x != ""])

  # Replace NAs or missing sequences with "-"
  gc.list <- lapply(gc.list, function(x) if (any(is.na(x))) "-" else x)

  # Check for "*" in the sequences and process to reverse operons
  te <- lapply(gc.list, function(x) grep("\\*", x, value = TRUE))
  if (length(te) > 0 && any(sapply(te, length) > 0)) {
    ye <- unlist(lapply(te, function(x) substr(x[1], 1, 1)))
    torev <- which(ye == "<")
  } else {
    torev <- NULL
  }

  # Process if torev exists
  if (!is.null(torev) && length(torev) > 0) {
    te <- gc.list[torev]
    te <- lapply(te, function(x) gsub("<-|->", "", x))
    te <- lapply(te, rev)

    # Split sequences with "=" for further processing
    witheq <- grep("=", te)
    withouteq <- setdiff(seq_along(te), witheq)

    # Process sequences with "=" using straightenOperonSeq
    if (length(witheq) > 0) {
      ge <- te[witheq]
      ge <- lapply(ge, straightenOperonSeq)
      te[witheq] <- ge
    }

    # Handle sequences without "=" by adding "->" suffix
    if (length(withouteq) > 0) {
      ye <- lapply(te[withouteq], function(x) paste(x, "->", sep = ""))
      te[withouteq] <- ye
    }

    gc.list[torev] <- te
  }

  # Reassemble GenContext from gc.list and update prot data frame
  rev.gencontext <- sapply(gc.list, function(x) paste(x, collapse = ""))
  rev.gencontext <- gsub("=", "|", rev.gencontext)
  prot$GenContext <- rev.gencontext

  return(prot)
}

# Helper function to straighten sequences based on equal signs
straightenOperonSeq <- function(w) {
  y <- rep(NA, length(w))
  d <- 1
  b <- grep("\\*", w)

  for (j in b:length(w)) {
    if (w[j] == "=") d <- d * (-1)
    y[j] <- ifelse(d == 1 && w[j] != "=", paste(w[j], "->", sep = ""),
      ifelse(d == -1 && w[j] != "=", paste("<-", w[j], sep = ""), "=")
    )
  }

  if (b > 1) {
    d <- 1
    for (j in (b - 1):1) {
      if (w[j] == "=") d <- d * (-1)
      y[j] <- ifelse(d == 1 && w[j] != "=", paste(w[j], "->", sep = ""),
        ifelse(d == -1 && w[j] != "=", paste("<-", w[j], sep = ""), "=")
      )
    }
  }
  return(y)
}

# Example usage:
prot <- data.frame(GenContext = c("A>B", "C<D", "E=F*G", "H>I"))
reversed_prot <- reverseOperonSeq(prot)
print(reversed_prot)
