# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(rlang))
# suppressPackageStartupMessages(library(igraph))
# suppressPackageStartupMessages(library(visNetwork))
# conflicted::conflict_prefer("filter", "dplyr")

###########################
## GC Undirected Network ##
###########################

#' Domain Network
#'
#' @description
#' This function creates a domain network from the 'DomArch' column.
#'
#' A network of domains is returned based on shared domain architectures.
#'
#'
#' @param prot A data frame that contains the column 'DomArch'.
#' @param column Name of column containing Domain architecture from which nodes and edges are generated.
#' @param domains_of_interest
#' @param cutoff_type Character. Used to determine how data should be filtered. Either
#' \itemize{\item "Lineage" to filter domains based off how many lineages the Domain architecture appears in
#' \item "Total Count" to filter off the total amount of times a domain architecture occurs }
#' @param cutoff Integer. Only use domains that occur at or above the cutoff for total counts if cutoff_type is "Total Count".
#' Only use domains that appear in cutoff or greater lineages if cutoff_type is Lineage.
#' @param layout Character. Layout type to be used for the network. Options are:
#' \itemize{\item "grid" \item "circle" \item "random" \item "auto"}
#'
#' @importFrom dplyr filter select
#' @importFrom grDevices rainbow
#' @importFrom igraph E graph_from_edgelist layout.auto layout.circle layout_on_grid layout_randomly  plot.igraph V
#' @importFrom stringr str_replace_all str_split
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' domain_network(pspa)
#' }
createUndirectedGenomicContextNetwork <- function(prot, column = "GenContext", domains_of_interest, cutoff_type = "Lineage", cutoff = 1, layout = "grid") {
    # by domain networks or all, as required.
    # ye is either all of prot.list or centered on one domain

    column_name <- sym(column)
    if (cutoff_type == "Lineage") {
        lin_summary <- prot %>%
            summarizeDomArch_ByLineage() %>%
            summarizeDomArch()
        doms_above_cutoff <- (lin_summary %>% filter(totallin >= cutoff))[[column]]
    } else if (cutoff_type == "Total Count") { # Change this type?
        GC_above_cutoff <- (prot %>% totalGenContextOrDomArchCounts(column = column, cutoff = cutoff))[[column]]
    }

    prot <- prot[which(prot[[as_string(column_name)]] %in% GC_above_cutoff), ]

    # prot.list=total$arch
    # prot.list <- unlist(lapply(prot.list, function(x) gsub(x,pattern = "\\?",replacement = "X")))
    prot$GC.ntwrk <- as_vector(prot %>% select({{ column }})) %>% # testing with prot$DomArch.orig
        str_replace_all(coll(pattern = "\\?", ignore_case = T), "X")

    ### Do I need a GC of interest?
    # ye <- prot %>%
    #   dplyr::filter(grepl(pattern=domains_of_interest,
    #                       x=DomArch.ntwrk,
    #                       ignore.case=T))

    ## Split the GC column by arrows and '||'
    GC_split <- prot$GC.ntwrk %>% str_split(pattern = "<-|->|\\|\\|")


    # GC.list <- GC_split
    if (any(lengths(GC_split) == 1)) {
        GC_split <- GC_split[-which(lengths(GC_split) == 1)]
    }
    lead <- unlist(lapply(GC_split, function(x) sapply(1:(length(x) - 1), function(y) x[y]))) # list elements 1 through n-1
    follow <- unlist(lapply(GC_split, function(x) sapply(1:(length(x) - 1), function(y) x[y + 1]))) # list elements 2 through n
    pwise <- cbind(lead, follow)
    if (any(pwise == "" | pwise == "?")) {
        pwise <- pwise[-which((pwise[, 1] == "") | (pwise[, 1] == "?") | (pwise[, 2] == "") | (pwise[, 2] == "?")), ]
    }
    pwise.ass <- sapply(1:length(pwise[, 1]), function(x) paste(pwise[x, ], collapse = "->"))
    e.sz <- sort(table(pwise.ass), decreasing = T) # Edge weights
    v.sz <- sort(table(pwise), decreasing = T) # Vertex weights
    pwise <- strsplit(names(e.sz), "\\->")
    pwise <- cbind(unlist(lapply(pwise, function(x) x[1])), unlist(lapply(pwise, function(x) x[2])))

    # Plotting the network
    g <- graph_from_edgelist(pwise, directed = TRUE)
    # scaling vertex size
    V(g)$size <- v.sz[V(g)$name]
    V(g)$size <- (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * 20 + 10 # scaled by degree
    # setting vertex color by size
    V(g)$color <- rainbow(5, alpha = .5)[round((V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * 4 + 1)]
    V(g)$frame.color <- V(g)$color
    # scaling edge width
    E(g)$width <- e.sz
    E(g)$width <- ifelse(log(E(g)$width) == 0, .3, log(E(g)$width))
    # coloring edges by width
    ew <- c(2.7, 4.5)
    E(g)$color <- sapply(
        E(g)$width,
        function(x) if (x >= ew[1] && x <= ew[2]) E(g)$color <- "cadetblue" else if (x > ew[2]) E(g)$color <- "maroon" else E(g)$color <- "gray55"
    )


    switch(layout,
        "random" = plot.igraph(g, layout = layout_randomly(g), vertex.label.dist = 0, asp = 0, edge.curved = F, vertex.label.color = "black"),
        "grid" = plot.igraph(g, layout = layout_on_grid(g), vertex.label.dist = 0, asp = 0, edge.curved = F, vertex.label.color = "black"),
        "circle" = plot.igraph(g, layout = layout.circle(g), vertex.label.dist = 0, asp = 0, edge.curved = F, vertex.label.color = "black"),
        "auto" = plot.igraph(g, layout = layout.auto(g), vertex.label.dist = 0, asp = 0, edge.curved = F, vertex.label.color = "black")
    )
}


#########################
## GC Directed Network ##
#########################
#' Genomic Context Directed Network
#'
#' @description
#' This function creates a Genomic Context network from the 'GenContext' column.
#'
#' A network of Genomic Context is returned.
#'
#'
#' @param prot A data frame that contains the column 'GenContext'.
#' @param domains_of_interest Character vector of domains of interest.
#' @param column Name of column containing Genomic Context from which nodes and edges are generated.
#' @param cutoff Integer. Only use GenContexts that occur at or above the cutoff percentage for total count
#' @param layout Character. Layout type to be used for the network. Options are:
#' \itemize{\item "grid" \item "circle" \item "random" \item "auto" \item "nice"}
#' @param directed Is the network directed?
#'
#' @importFrom dplyr distinct filter group_by pull select summarize
#' @importFrom purrr as_vector map map2
#' @importFrom rlang sym
#' @importFrom stringr str_replace_all
#' @importFrom visNetwork visIgraphLayout visLegend visNetwork visOptions
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gc_directed_network(pspa, column = "GenContex", cutoff = 55)
#' }
createGenomicContextNetwork <- function(prot, domains_of_interest, column = "GenContext",
    cutoff = 40,
    layout = "grid",
    directed = TRUE) {
    column_name <- sym(column)


    # Perform cutoff through totalGenContextOrDomArchCounts
    prot_tc <- prot %>% totalGenContextOrDomArchCounts(column = column, cutoff = cutoff)

    within_list <- prot_tc %>%
        select({{ column_name }}) %>%
        distinct()
    within_list <- pull(within_list, {{ column_name }})

    prot <- prot %>% filter({{ column_name }} %in% within_list)

    # Below is obsolete kinda
    prot$DomArch.ntwrk <- as_vector(prot %>% select({{ column }})) %>%
        str_replace_all(coll(pattern = "\\?", ignore_case = T), "X")

    # # Make sure all observations contain the domain of interest
    # domains_of_interest_regex = paste(domains_of_interest, collapse = "|")
    # ye <- prot %>%
    #   dplyr::filter(grepl(pattern=domains_of_interest_regex,
    #                       x=DomArch.ntwrk,
    #                       ignore.case=T))

    # ##Separating column and converting to atomic vector prevents coercion
    # ye <- ye$DomArch.ntwrk  %>% str_split(pattern="\\+")
    #
    # GC.list <- ye
    # if(any(lengths(GC.list)==1)) {
    #   GC.list=GC.list[-which(lengths(GC.list)==1)]
    # }


    #### rev_gc(gc_row)
    #### Reverse the order of the split up GC elements whenever a change in transcription direction is found?
    rev_gc <- function(x) {
        direction <- 0
        op_dir <- "<-|->"
        # if direction is 1, '->' are the correct arrows
        # if direction is -1, '<-'
        last <- ""
        for (i in 1:length(x)) {
            if (direction == 0 && grepl("<-|->", x[i])) {
                # Set initial direction of the GC
                if (grepl("^<-", x[i])) {
                    op_dir <- "->"
                    direction <- -1
                } else {
                    op_dir <- "<-"
                    direction <- 1
                }
            } else if (grepl(op_dir, x[i])) {
                # Arrow going in the opposite direction has been found
                # Reverse order of the remaining list--- Do I need to add a " " before this?
                x[i:length(x)] <- rev(x[i:length(x)])
                if (direction == 1) {
                    op_dir <- "->"
                    direction <- -1
                } else if (direction == -1) {
                    op_dir <- "<-"
                    direction <- 1
                }
            }
        }
        return(x)
    }

    # Get domain counts before eliminating domarchs with no edges
    wc <- elements2Words(prot = prot, column = column, conversion_type = "gc2da") %>% words2WordCounts()
    nodes <- data.frame(id = wc$words, label = wc$words, size = wc$freq)

    max_size <- max(nodes$size)
    min_size <- min(nodes$size)
    nodes <- nodes %>% mutate(size = (size - min_size) / ((max_size - min_size)) * 20 + 10)
    max_font_size <- 43
    nodes <- nodes %>% mutate(font.size = purrr::map(size, function(x) min(x * 2, max_font_size)))

    max_size <- max(nodes$size)
    min_size <- min(nodes$size)
    nodes <- nodes %>% mutate(color = unlist(purrr::map(
        nodes$size,
        function(x) {
            return(rainbow(5, alpha = .5)[round((x - min_size) / (max_size - min_size) * 4 + 1)])
        }
    )))
    nodes$frame.color <- nodes$color



    gencontext <- prot$DomArch.ntwrk

    # Split GC
    gencontext <- gsub(pattern = ">", replacement = ">|", x = gencontext)
    gencontext <- gsub(pattern = "<", replacement = "|<", x = gencontext)

    gen_split <- gencontext %>% strsplit("\\|")

    # Remove non edges
    if (any(lengths(gen_split) == 1)) {
        gen_split <- gen_split[-which(lengths(gen_split) == 1)]
    }

    # Apply rev_gc to all rows
    gen_split_reoriented <- map(.x = gen_split, .f = rev_gc)

    # Remove all arrows
    gen_split_reoriented <- map2(.x = gen_split_reoriented, .y = "<-|->", .f = str_remove_all)

    from <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x) - 1), function(y) x[y]))) # list elements 1 through n-1
    to <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x) - 1), function(y) x[y + 1]))) # list elements 2 through n

    # Edge list
    pwise <- cbind(from, to)

    # Remove rows that have bad characters
    if (any(pwise == "" | pwise == "X" | pwise == "X(s)")) {
        pwise <- pwise[-which((pwise[, 1] == "") | (pwise[, 1] == "X") | (pwise[, 1] == "X(s)") | (pwise[, 2] == "") | (pwise[, 2] == "X") | (pwise[, 2] == "X(s)")), ]
    }

    edges <- data.frame(from = pwise[, 1], to = pwise[, 2]) %>%
        group_by(from, to) %>%
        summarize(width = n())
    edges <- edges %>% mutate(width = ifelse(width == 1, .3, log(width)))
    ew <- c(2.7, 4.5)

    ColorEdges <- function(x) {
        if (x >= ew[1] && x <= ew[2]) {
            adjustcolor("cadetblue", alpha.f = .7)
        } else if (x > ew[2]) {
            adjustcolor("maroon", alpha.f = .5)
        } else {
            "gray55"
        }
    }

    edges <- edges %>% mutate(color = unlist(purrr::map(width, ColorEdges)))

    if (directed) {
        vg <- visNetwork(nodes, edges, width = "100%", height = "600px") %>%
            visEdges(arrows = "to", smooth = T)
    } else {
        vg <- visNetwork(nodes, edges, width = "100%", height = "600px") %>%
            visEdges(smooth = T)
    }

    vg <- vg %>%
        visOptions(highlightNearest = TRUE) %>%
        visLegend()

    vg <- switch(layout,
        "nice" = visIgraphLayout(vg, "layout_nicely"),
        "random" = visIgraphLayout(vg, "layout_randomly"),
        "grid" = visIgraphLayout(vg, "layout_on_grid"),
        "circle" = visIgraphLayout(vg, "layout.circle"),
        "auto" = visIgraphLayout(vg, "layout.auto")
    )
    vg
}
