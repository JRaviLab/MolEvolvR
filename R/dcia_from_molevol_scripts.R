#' Plot Protein Length Distribution
#'
#' This function plots the distribution of protein lengths from the provided data.
#' 
#' @references 
#' \url{https://journals.asm.org/doi/epub/10.1128/jb.00163-22}
#' Figure 2A
#'
#' @param data A data frame containing protein information.
#' @param query_data A data frame containing protein query information.
#' @param superkingdom_filter A character vector specifying the superkingdom(s) to filter the data by. Valid
#' entries include "bacteria", "virus", "eukaryota", "archaea" or a combination. 
#' Default is NULL, meaning no filtering.
#' @param query_filter A character vector specifying one or many Query Proteins.
#' Default is NULL, meaning no filtering.
#' @param query_lineage_filter A character vector specifying one or many Query Lineage Proteins.
#' Default is NULL, meaning no filtering. 
#' @param show.median Logical. Whether or not to draw a median line on the plot
#' @param show.percentiles Logical. Whether or not to draw 25th and 75th percentile lines on the plot.
#'
#' @return A ggplot object showing the distribution of protein lengths.
#' @export
#'
#' @examples
#' \dontrun{
#' plotProteinLengthDist(data, superkingdom_filter = "bacteria")
#' }

# Plot protein length distribution as a box plot
plotProteinLengthDist_byLineage <- function(data, query_data, superkingdom_filter = NULL, query_filter = NULL, query_lineage_filter = NULL, show.median = TRUE, show.percentiles = TRUE) {
  ## Verify Required Subject Length column (`SLength`) is present in data
  slength_col_check <- 'SLength' %in% colnames(data)
  if(slength_col_check == FALSE) {
    rlang::abort('SLength column is not present in supplied data and is required for this plot.', class = 'data-error')
    } 
  ## Create Regex Filters Based on user input
  ### Superkingdom
  if(any(!str_detect(superkingdom_filter, regex('bacteria|virus|eukaryota|archaea', ignore_case = TRUE)))) {
    rlang::abort(message = 'Please filter by "bacteria", "virus", "eukaryota", "archaea" or a combination as a vector')
    } else {
      superkingdom_filter_regex <- str_trim(glue::glue_collapse(superkingdom_filter, sep = '|'))
      }
  ### Query Protein
  query_filter_regex <- str_trim(glue::glue_collapse(query_filter, sep = '|'))
  query_lineage_filter_regex <- str_trim(glue::glue_collapse(query_lineage_filter, sep = '|'))

  ### Query Data Prep
  query_data_filter <- 
    query_data %>% 
    distinct(QueryName, Lineage)

  
  ## Data Prep: 
  ## - create superkingdom column, parsing from Lineage
  ## - remove empty or NA Lineage values
  df <- 
    data %>% 
    mutate(superkingdom = str_extract(Lineage, regex('^[^>]+'))) %>% 
    filter(Lineage != "") %>% 
    drop_na(Lineage)

  ## Apply Filters: Always remove 'candidatus|candidate' along with any user specified filtering
  df_filtered <- 
    df %>%
    filter(!str_detect(Lineage, regex('candidatus|candidate', ignore_case = TRUE))) %>% 
    { if(is.null(superkingdom_filter)) . else filter(., str_detect(Lineage, regex(superkingdom_filter_regex, ignore_case = TRUE))) } %>% 
    { if(is.null(query_filter)) . else filter(., str_detect(QueryName, regex(query_filter_regex, ignore_case = TRUE))) } %>% 
    distinct(AccNum, .keep_all = TRUE) 
  
  df_filtered_long <- 
    df_filtered %>%
    pivot_longer(c(SLength, SStart, SEnd, QLength, QStart, QEnd), names_to = c("type", 'attribute'), names_sep = 1, values_to = "length") %>% 
    mutate(attribute = fct_relevel(attribute, c("Length", "Start", "End"))) %>% 
    filter(attribute == 'Length')

  ## Calculate Stats: Median, 25th and 75th percentiles
  df_filtered_stats <- 
    df_filtered_long %>% 
    {if(is.null(query_lineage_filter)) . else filter(., str_detect(Lineage, regex(query_lineage_filter_regex, ignore_case = TRUE))) } %>% 
    group_by(type, attribute) %>%
    mutate(median_length = median(length), 
           first_quartile_length = quantile(length, c(0.25)), 
           third_quartile_length = quantile(length, c(0.75))
          )
  ## Create Plot
  if(nrow(df_filtered_long) > 0){
    ProteinLengthDist <- 
      df_filtered_long %>% 
      filter(type == 'S') %>% 
      ggplot(aes(x=Lineage, y = length, fill = superkingdom)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.color = "black") +
      scale_fill_viridis(discrete = TRUE, alpha=0.6, guide = 'none') +
      labs(y = "Length (in amino acids)", x = "Lineage(s)") +
      theme_minimal() +
      theme(text = element_text(size = 18),
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.8, 
                                       hjust=0
                                      ),
            legend.title = element_text(family = "Inter", face = "bold"),
            legend.text = element_text(family = "Inter"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect("white", linewidth = 0),
            panel.border = element_blank()) +
      geom_point(data = df_filtered_stats %>% filter(type == 'Q') %>% inner_join(query_data_filter) %>% distinct(), aes(color = QueryName, shape = QueryName), size = 3) +
      scale_color_brewer(palette = 'Dark2') +
      scale_x_discrete(position = "top") +
      guides(color = guide_legend("Query Protein"), shape = guide_legend("Query Protein"))
    if(show.median & show.percentiles){
      ProteinLengthDist +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = median_length), color = "blue") +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = first_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5) +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = third_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5)
      } else if (show.median) {
        ProteinLengthDist +
          geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = median_length), color = "blue")
        } else if (show.percentiles){
          ProteinLengthDist +
            geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = first_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5) +
            geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = third_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5)
          } else {
            ProteinLengthDist
            }
    } else {
      rlang::abort(message = 'No data returned with configured Superkingdom/Query filter(s)', class = 'filter-error')
      }
}

# Plot protein length distribution as a box plot
plotProteinLengthDist_byQuery <- function(data, query_data, superkingdom_filter = NULL, query_filter = NULL, query_lineage_filter = NULL, show.median = TRUE, show.percentiles = TRUE) {
  ## Verify Required Subject Length column (`SLength`) is present in data
  slength_col_check <- 'SLength' %in% colnames(data)
  if(slength_col_check == FALSE) {
    rlang::abort('SLength column is not present in supplied data and is required for this plot.', class = 'data-error')
    } 
  ## Create Regex Filters Based on user input
  ### Superkingdom
  if(any(!str_detect(superkingdom_filter, regex('bacteria|virus|eukaryota|archaea', ignore_case = TRUE)))) {
    rlang::abort(message = 'Please filter by "bacteria", "virus", "eukaryota", "archaea" or a combination as a vector')
    } else {
      superkingdom_filter_regex <- str_trim(glue::glue_collapse(superkingdom_filter, sep = '|'))
      }
  ### Query Protein
  query_filter_regex <- str_trim(glue::glue_collapse(query_filter, sep = '|'))
  query_lineage_filter_regex <- str_trim(glue::glue_collapse(query_lineage_filter, sep = '|'))

  ### Query Data Prep
  query_data_filter <- 
    query_data %>% 
    distinct(QueryName, Lineage)

  
  ## Data Prep: 
  ## - create superkingdom column, parsing from Lineage
  ## - remove empty or NA Lineage values
  df <- 
    data %>% 
    mutate(superkingdom = str_extract(Lineage, regex('^[^>]+'))) %>% 
    filter(Lineage != "") %>% 
    drop_na(Lineage)

  ## Apply Filters: Always remove 'candidatus|candidate' along with any user specified filtering
  df_filtered <- 
    df %>%
    filter(!str_detect(Lineage, regex('candidatus|candidate', ignore_case = TRUE))) %>% 
    { if(is.null(superkingdom_filter)) . else filter(., str_detect(Lineage, regex(superkingdom_filter_regex, ignore_case = TRUE))) } %>% 
    { if(is.null(query_filter)) . else filter(., str_detect(QueryName, regex(query_filter_regex, ignore_case = TRUE))) } %>% 
    distinct(AccNum, .keep_all = TRUE) 
  
  df_filtered_long <- 
    df_filtered %>%
    pivot_longer(c(SLength, SStart, SEnd, QLength, QStart, QEnd), names_to = c("type", 'attribute'), names_sep = 1, values_to = "length") %>% 
    mutate(attribute = fct_relevel(attribute, c("Length", "Start", "End"))) %>% 
    filter(attribute == 'Length')

  ## Calculate Stats: Median, 25th and 75th percentiles
  df_filtered_stats <- 
    df_filtered_long %>% 
    {if(is.null(query_lineage_filter)) . else filter(., str_detect(Lineage, regex(query_lineage_filter_regex, ignore_case = TRUE))) } %>% 
    group_by(type, attribute) %>%
    mutate(median_length = median(length), 
           first_quartile_length = quantile(length, c(0.25)), 
           third_quartile_length = quantile(length, c(0.75))
          )
  ## Create Plot
  if(nrow(df_filtered_long) > 0){
    ProteinLengthDist <- 
      df_filtered_long %>% 
      filter(type == 'S') %>% 
      ggplot(aes(x=Query, y = length, fill = superkingdom)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.color = "black") +
      scale_fill_viridis(discrete = TRUE, alpha=0.6, guide = 'none') +
      labs(y = "Length (in amino acids)", x = "Query") +
      theme_minimal() +
      theme(text = element_text(size = 18),
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.8, 
                                       hjust=0
                                      ),
            legend.title = element_text(family = "Inter", face = "bold"),
            legend.text = element_text(family = "Inter"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect("white", linewidth = 0),
            panel.border = element_blank()) +
      # geom_point(data = df_filtered_stats %>% filter(type == 'Q') %>% inner_join(query_data_filter) %>% distinct(), aes(color = QueryName, shape = QueryName), size = 3) +
      scale_color_brewer(palette = 'Dark2') +
      scale_x_discrete(position = "top") +
      guides(color = guide_legend("Query Protein"), shape = guide_legend("Query Protein"))
    if(show.median & show.percentiles){
      ProteinLengthDist +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = median_length), color = "blue") +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = first_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5) +
        geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = third_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5)
      } else if (show.median) {
        ProteinLengthDist +
          geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = median_length), color = "blue")
        } else if (show.percentiles){
          ProteinLengthDist +
            geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = first_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5) +
            geom_hline(data = df_filtered_stats %>% filter(type == 'S'), aes(yintercept = third_quartile_length), linetype = 'dashed', color = "light blue", alpha = 0.5)
          } else {
            ProteinLengthDist
            }
    } else {
      rlang::abort(message = 'No data returned with configured Superkingdom/Query filter(s)', class = 'filter-error')
      }
}