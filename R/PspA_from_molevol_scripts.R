label_creator <- function(label_text, font_size) {
  ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y, label = label_text)) +
  geom_text(size = font_size, fontface = "bold", family = "Inter") +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.title = element_text(hjust = 0.5)
  )
}

#' @title Stacked Lineage Plot by Domain Architecture
#' @description Renders a stacked bar plot of lineage occurrences of the top domain
#' architectures (shown on the left); similar to the [PSP paper, Figure 4](https://doi.org/10.1128/msystems.00847-23)
#' @param prot Data frame containing DomArch (domain architecture; default: Pfam) and Lineage columns
#' @param column Column that contains domain architectures in the prot dataframe (default: Pfam)
#' @param cutoff Minimum fraction of lineage-wise occurrences
#' @param Lineage_col Column containing the lineages in the prot dataframe
#' @param reduce_lineage Abbreviate lineages to a shorter number of characters
#' @param label.size Size of the non-legend text labels, e.g. the lineage tick marks on the x axis
#' @param legend.position Position of the legend
#' @param legend.text.size Size of the legend text
#' @param legend.cols Number of columns in the legend
#' @param legend.size Size of the legend key
#' @param char_width Width of each DomArch character to size its enclosing box appropriately
#' @param coord_flip Whether to flip the plot so the x-axis is number of lineage occurrences (TRUE/FALSE)
#' @param legend Whether to show the legend (TRUE/FALSE)
#' @param verbose Whether to print verbose output, e.g. the number of rows in the chart, unique colors, etc. (TRUE/FALSE)
#' @param base_axes_font_size Font size for the labels at the bottom of the chart
#' @return A list with the following keys: 'plot', the ggplot2 object, and
#' 'chart_rows', the number of rows in the chart
stacked_lin_plot <- function(prot, column = "DomArch", cutoff, Lineage_col = "Lineage",
                             xlabel = "Domain Architecture",
                             reduce_lineage = TRUE,
                             label.size = 8,
                             legend.position = c(0.7, 0.4),
                             legend.text.size = 10,
                             legend.cols = 2,
                             legend.size = 0.7,
                             char_width = 0.2,
                             coord_flip = TRUE, legend = TRUE, verbose = FALSE,
                             base_axes_font_size = 6) {
  CPCOLS <- c(
    "#AFEEEE", "#DDA0DD", "#EE2C2C", "#CDBE70", "#B0B099",
    "#8B2323", "#EE7600", "#EEC900", "chartreuse3", "#0000FF",
    "#FFD900", "#32CD32", "maroon4", "cornflowerblue", "darkslateblue",
    "#AB82FF", "#CD6889", "#FFA07A", "#FFFF00", "#228B22",
    "#FFFFE0", "#FFEC8B", "peru", "#668B8B", "honeydew",
    "#A020F0", "grey", "#8B4513", "#191970", "#00FF7F",
    "lemonchiffon", "#66CDAA", "#5F9EA0", "#A2CD5A", "#556B2F",
    "#EEAEEE", "thistle4", "#473C8B", "#FFB6C1", "#8B1C62",
    "#FFE4B5", "black", "#FF7F50", "#FFB90F", "#FF69B4", "#836FFF",
    "#757575", "#CD3333", "#EE7600", "#CDAD00", "#556B2F", "#7AC5CD"
  )
  col <- sym(column)

  if (reduce_lineage) {
    prot <- shorten_lineage(prot, Lineage_col, abr_len = 3)
  }

  total_count <- total_counts(prot, column, cutoff, lineage_col = Lineage_col)

  # Order bars by descending freq
  order <- total_count %>%
    select({{ col }}) %>%
    unique() %>%
    pull({{ col }})
  if (coord_flip) {
    order <- order %>% rev()
  }

  prot_data <- total_count
  prot_data[[Lineage_col]] <- unlist(map(prot_data[[Lineage_col]], function(x) {
    gt_pos <- gregexpr(pattern = ">", x)[[1]][2]
    # If there is not a second '>' in the lineage
    if (is.na(gt_pos)) {
      x
    } else {
      x
      # substr(x,1, (gt_pos-1) )
    }
  }))

  prot_data[, column] <- factor(pull(prot_data, {{ col }}))

  # sort the bar chart descending by the total number of lineage occurrences
  prot_data[[column]] <- reorder(prot_data[[column]], prot_data$totalcount)

  # replace any number of NAs with a single NA
  prot_data[[Lineage_col]] = prot_data[[Lineage_col]] %>%
   str_replace(pattern = "^(NA)+$", replacement = "NA")

  # ---------------------------------------------
  # start creating graphical elements
  # ---------------------------------------------

  # the full graphic consists of the following 4 panes:
  # 1a. the list of domains on the left
  # 2a. the stacked vertical bar chart on the right
  # 1b. a "Domain Architecture" label on the bottom left
  # 2b. a "Number of Occurrences" label on the bottom right

  if (coord_flip) {
    if (legend) {
      stacked_bar <- ggplot(prot_data, aes_string(fill = Lineage_col, y = "count", x = column)) +
        geom_bar(position = "stack", stat = "identity", color = "white", width=1) +
        coord_flip() +
        xlab("Group") +
        ylab("Number of proteins") +
        theme_minimal() +
        scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
        scale_x_discrete(expand = c(0, 0)) +  # Remove space on the x-axis
        scale_y_continuous(expand = c(0, 0)) +  # Remove space on the y-axis
        theme(
          legend.position = legend.position,
          # legend.background = element_rect(fill = "white", color = "white"),
          legend.background = element_blank(),
          legend.text = element_text(size = legend.text.size, family = "Inter"),
          legend.title = element_text(size = legend.text.size + 2, family = "Inter"),
          legend.key.size = unit(legend.size, "cm"),
          # legend.key.height = unit(2, "cm"),
          # legend.key.width = unit(0.9, "cm"),
          legend.spacing = unit(0.4, "cm"),
          legend.spacing.x = unit(1.5, "cm"), # makes lineages not overlap
          axis.text = element_text(size = label.size, angle = 315, hjust = 1, vjust = 0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(), # element_line(colour = "black"),
          plot.background = element_rect("white", linewidth = 0),
          panel.border = element_blank(),
          # axis.text.x = ele

          axis.text.x = element_text(hjust=0, face="bold", family = "Inter"),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),

          plot.margin = margin(0, 0, 0, 0)
        ) +
        guides(fill = guide_legend(ncol = legend.cols, reverse = TRUE))
    } else {
      stacked_bar <- ggplot(prot_data, aes_string(fill = Lineage_col, y = "count", x = column)) +
        geom_bar(position = "stack", stat = "identity", color = "white") +
        coord_flip() +
        xlab("Group") +
        ylab("Number of proteins") +
        theme_minimal() +
        scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
        theme(
          legend.position = "none",
          legend.spacing = unit(0.4, "cm"),
          axis.text = element_text(size = label.size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.background = element_rect("white", linewidth = 0),
          panel.border = element_blank()
          # axis.text.x = ele
        )
    }
  } else {
    stacked_bar <- ggplot(prot_data, aes(fill = {{ Lineage_col }}, y = count, x = {{ col }})) +
      geom_bar(position = "stack", stat = "identity") +
      xlab("") +
      ylab("") +
      theme_minimal() +
      theme(
        legend.position = legend.position,
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = legend.text.size),
        legend.title = element_text(size = legend.text.size + 2),
        legend.key.size = unit(legend.size, "cm"),
        axis.text = element_text(size = label.size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        # axis.text.x = ele
      ) +
      guides(fill = guide_legend(
        ncol = legend.cols, reverse = TRUE
      ))
  }

  # ---
  # make the list of domain architectures
  # ---

  prot_data$labels <- reorder(prot_data[[column]], prot_data$totalcount)

  # Step 1: Split labels on '+'
  label_components <- prot_data %>%
    distinct(labels) %>%
    mutate(id = row_number()) %>%  # Add an identifier for each row
    separate_rows(labels, sep = "\\+") %>%  # Split on '+' and create new rows
    group_by(id) %>%
    mutate(
      labels = rev(labels),
      width = nchar(labels) * char_width,
      total_width = sum(width),
      position = cumsum(dplyr::lag(width, default=0)),
      boxmax = position + width
    ) %>%
    ungroup()

  # get the number of rows in the resulting chart
  chart_rows <- nrow(label_components %>% distinct(id))

  # figure out how many values we need and generate the colors
  # (this is to handle cases where there are way more domains than there
  # are colors in the typical palettes)
  unique_colors <- length(unique(label_components$labels))
  hcl_palette <- colorspace::qualitative_hcl(
    n = unique_colors, palette = "Set3"
  )

  # Step 2: Text Labels with Multiple Boxes
  text_labels <- ggplot(label_components, aes(
      x = position, y = id, label = labels, fill = labels,
      xmin = position, xmax = boxmax,
      ymin = id, ymax = id + 1
    )) +
    geom_rect(show.legend = FALSE, fill = "white") +
    geom_rect(show.legend = FALSE, color = "#525252", aes(
      xmin = position, xmax = boxmax,
      ymin = id + 0.2, ymax = id + 1 - 0.2
    )) +
    geom_fit_text(
      family = "Inter",
      aes(
        xmin = position, xmax = boxmax,
        ymin = id + 0.2, ymax = id + 1 - 0.2,
        label = labels
      ),
      min.size = 2, grow = FALSE, show.legend = FALSE
    ) +
    scale_y_reverse(expand=c(0,0)) + # make domain sets read in the correct order, remove all y-axis padding 
    scale_x_reverse(expand=c(0.0125, 0)) + # make domain sets right-justified, with some padding
    scale_fill_manual(values = hcl_palette, na.value = "#A9A9A9", aesthetics = c("colour", "fill")) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, label.size + 21, 0)
    )

  # Step 3. Assembly of the final plot with multiple subplots

  # create text labels under each figure
  left_title <- label_creator("Domain Architectures", font_size = base_axes_font_size)
  right_title <- label_creator("Number of Occurrences per Lineage", font_size = base_axes_font_size)

  # combine the plots side by side, with their text labels below
  fullplot <- grid.arrange(
    text_labels, stacked_bar, left_title, right_title,
    ncol = 2, widths = c(3, 5), heights = c(5, unit(0.4, "inches"))
  )

  return(list(
    "plot" = fullplot,
    "chart_rows" = chart_rows
  ))
}
