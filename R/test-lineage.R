# Load necessary packages
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

# Global LevelReduction function
LevelReduction <- function(lin, level = 1) {
    if (level == 1) {
        gt_loc <- str_locate(lin, ">")[[1]]
        if (is.na(gt_loc)) {
            return(lin)
        } else {
            lin <- substring(lin, first = 0, last = (gt_loc - 1))
            return(lin)
        }
    }
    gt_loc <- str_locate_all(lin, ">")[[1]]
    l <- length(gt_loc) / 2
    if (level > l) {
        return(lin)
    } else {
        gt_loc <- as.numeric(gt_loc[level, 1])
        lin <- substring(lin, first = 0, last = (gt_loc - 1))
        return(lin)
    }
}

# Global GetKingdom function
GetKingdom <- function(lin) {
    gt_loc <- str_locate(lin, ">")[ "start"]
    if (is.na(gt_loc)) {
        return(lin)
    } else {
        kingdom <- substring(lin, first = 0, last = (gt_loc - 1))
        return(kingdom)
    }
}

# Refactored LineagePlot function
LineagePlot <- function(prot, domains_of_interest, level = 3, label.size = 8) {
    all_grouped <- data.frame("Query" = character(0), "Lineage" = character(0), "count" = integer())

    for (dom in domains_of_interest) {
        domSub <- prot %>% filter(grepl(dom, GenContext, ignore.case = TRUE))
        domSub <- domSub %>%
            group_by(Lineage) %>%
            summarize("count" = n(), .groups = "drop")  # Add .groups argument

        domSub$Query <- dom
        all_grouped <- dplyr::union(all_grouped, domSub)
    }

    all_grouped <- all_grouped %>% mutate(ReducedLin = unlist(purrr::map(Lineage, LevelReduction)))

    all_grouped_reduced <- all_grouped %>%
        group_by(Query, ReducedLin) %>%
        summarize("count" = sum(count), .groups = "drop") %>%  # Add .groups argument
        mutate(Kingdom = unlist(purrr::map(ReducedLin, GetKingdom)))

    lin_counts <- all_grouped_reduced %>%
        group_by(Kingdom, ReducedLin) %>%
        summarize("count" = n(), .groups = "drop")  # Add .groups argument

    bacteria_colors <- rep("#d94f25", length(grep("bacteria", lin_counts$Kingdom)))
    archaea_colors <- rep("#26662d", length(grep("archaea", lin_counts$Kingdom)))
    eukaryota_colors <- rep("#0000ff", length(grep("eukaryota", lin_counts$Kingdom)))
    virus_colors <- rep("#4f4f4f", length(grep("virus", lin_counts$Kingdom)))

    colors <- append(archaea_colors, bacteria_colors) %>%
        append(eukaryota_colors) %>%
        append(virus_colors)

    all_grouped_reduced$ReducedLin <- map(
        all_grouped_reduced$ReducedLin,
        function(lin) {
            gt_loc <- str_locate(lin, ">")[ "start"]
            if (is.na(gt_loc)) {
                return(lin)
            } else {
                lin <- substring(lin, first = (gt_loc + 1), last = 100)
                return(lin)
            }
        }
    ) %>% unlist()

    ordered_lin <- all_grouped_reduced %>% arrange(Kingdom)

    all_grouped_reduced$Query <- factor(
        x = all_grouped_reduced$Query,
        levels = rev(domains_of_interest)
    )
    all_grouped_reduced$ReducedLin <- factor(
        x = all_grouped_reduced$ReducedLin,
        levels = unique(ordered_lin$ReducedLin)
    )

    # ggplot for the final plot
    ggplot(
        data = all_grouped_reduced,
        aes(x = ReducedLin, y = Query)  # Use aes() instead of aes_string()
    ) +
        geom_tile(
            data = subset(all_grouped_reduced, !is.na(count)),
            aes(fill = count),
            colour = "darkred", linewidth = 0.3  # Use linewidth instead of size
        ) +
        scale_fill_gradient(low = "white", high = "darkred") +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_text(
                angle = 90, hjust = 1, vjust = 0.2,
                color = colors, size = label.size + 5
            ),
            axis.text = element_text(size = label.size)
        )
}

# Test code
prot <- data.frame(
    Lineage = c("Animalia>Chordata>Mammalia", "Plantae>Tracheophyta", "Fungi>Ascomycota", "Bacteria>Proteobacteria>Gammaproteobacteria"),
    GenContext = c("PspA", "PspB", "PspC", "PspD")
)

domains_of_interest <- c("PspA", "PspB", "PspC", "PspD")

# Call the plotting function
LineagePlot(prot, domains_of_interest, level = 2, label.size = 10)
