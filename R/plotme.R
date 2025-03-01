# Taken from https://github.com/yogevherz/plotme/blob/master/R/count_to_sunburst_treemap.R
#' Create an interactive plotly from count data
#'
#' @description
#' These functions help you quickly create interactive hierarchical plots
#' from categorical data. They expect the summary of the data created by
#' `dplyr::count()` and produce either a sunburst plot (`plotSunburst()`) or
#' a treemap plot (`plotTreemap()`)
#'
#' @param count_data An output of dplyr::count(), tibble or data frame
#' @param fill_by_n If TRUE, uses a continuous scale to fill plot by group size
#' @param sort_by_n If TRUE, sorts groups in plot by size, if FALSE sorts them alphabetically
#'
#' @importFrom plotly plot_ly
#' @importFrom purrr exec
#'
#' @export
#' @examples
#' library(dplyr)
#' starwars_count <- count(starwars, species, eye_color, name)
#'
#' # sunburst plot
#' plotSunburst(starwars_count)
#'
#' # fill by group size
#' plotSunburst(starwars_count, fill_by_n = TRUE)
#'
#' # treemap plot, ordered by group size
#' plotTreemap(starwars_count, sort_by_n = TRUE)
#'
#' # display al charchaters by homeworld
#' starwars %>%
#'     count(homeworld, name) %>%
#'     plotTreemap(sort_by_n = TRUE)
#'
plotSunburst <- function(count_data, fill_by_n = FALSE, sort_by_n = FALSE, maxdepth = 2) {
    params <- prepareColumnParams(count_data, fill_by_n, sort_by_n)

    purrr::exec(plotly::plot_ly,
        !!!params,
        type = "sunburst",
        branchvalues = "total"
    )
}


#' @param count_data A data frame containing the data.
#' @param fill_by_n Logical indicating if fill color is based on counts.
#' @param sort_by_n Logical indicating if data should be sorted by counts.
#'
#' @importFrom plotly plot_ly
#' @importFrom purrr exec
#'
#' @export
#' @rdname plotSunburst
plotTreemap <- function(count_data, fill_by_n = FALSE, sort_by_n = FALSE) {
    params <- prepareColumnParams(count_data, fill_by_n, sort_by_n)

    purrr::exec(plotly::plot_ly,
        !!!params,
        type = "treemap",
        branchvalues = "total",
        hoverinfo = "text"
    )
}


#' prepareColumnParams
#'
#' @param count_data A data frame containing the data.
#' @param fill_by_n Logical indicating if fill color is based on counts.
#' @param sort_by_n Logical indicating if data should be sorted by counts.
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_rows mutate
#' @importFrom purrr map
#'
#' @return A data frame of parameters for treemap visualization.
#' @export
#'
#' @examples
#' \dontrun{
#' count_data <- data.frame(Category = c("A", "B", "C"),
#'                           n = c(10, 20, 15))
#' params <- prepareColumnParams(count_data, fill_by_n = TRUE, sort_by_n = FALSE)
#' params
#' }
prepareColumnParams <- function(count_data, fill_by_n, sort_by_n) {
    validateCountDF(count_data)
    assertthat::assert_that(is.logical(fill_by_n),
        length(fill_by_n) == 1,
        msg = "fill_by_n must be either TRUE or FALSE"
    )
    assertthat::assert_that(is.logical(sort_by_n),
        length(sort_by_n) == 1,
        msg = "sort_by_n must be either TRUE or FALSE"
    )

    count_data <- .all_non_n_cols_to_char(count_data)

    category_num <- ncol(count_data) - 1

    params <- purrr::map(1:category_num,
        prepareSingleColumnParams,
        df = count_data,
        root = ""
    ) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(sort = sort_by_n)

    if (fill_by_n) {
        params <- params %>%
            dplyr::mutate(marker = list(
                colorbar = list(
                    bgcolor = ""
                )
            ))
    }
    params
}

#' prepareSingleColumnParams
#'
#' @param df A data frame containing the data to be processed.
#' @param col_num An integer representing the column number to process.
#' @param root A string representing the root node for the treemap.
#'
#' @importFrom dplyr c_across group_by mutate rowwise select summarise ungroup
#' @importFrom stringr str_glue
#'
#' @return A data frame containing parameters for the specified column for
#' treemap visualization.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(Category = c("A", "A", "B", "B", "C"),
#'                  n = c(10, 20, 30, 40, 50))
#' params <- prepareSingleColumnParams(df, col_num = 1, root = "Root")
#' params
#' }
prepareSingleColumnParams <- function(df,
    col_num,
    root) {
    col_name <- names(df)[col_num]

    df %>%
        dplyr::group_by(dplyr::across(1:dplyr::all_of(col_num))) %>%
        dplyr::summarise(values = sum(.data$n), .groups = "drop") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            ids = paste(dplyr::c_across(1:!!col_num),
                collapse = ".->."
            ),
            parents = ifelse(!!col_num > 1,
                paste(dplyr::c_across(1:(!!col_num - 1)),
                    collapse = ".->."
                ),
                root
            )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            labels = .[[!!col_num]],
            hovertext = stringr::str_glue(
                "column: {col_name}\nvalue: {labels}\nn: {values}"
            )
        ) %>%
        dplyr::select(ids, parents, labels, values, hovertext)
}
#' validateCountDF
#'
#' @param var A data frame whose columns are to be converted.
#'
#' @importFrom assertthat assert_that has_name
#' @importFrom dplyr across mutate
#'
#' @return A data frame with non-'n' columns converted to character type.
#' @export
#'
#' @examples
#' \dontrun{
#' new_df <- .all_non_n_cols_to_char(my_data)
#' }
validateCountDF <- function(var) {
    msg <- paste(substitute(var), "must be a count dataframe (output of dplyr::count)")
    assertthat::assert_that(is.data.frame(var),
        assertthat::has_name(var, "n"),
        msg = msg
    )

    n_col <- var$n
    assertthat::assert_that(is.numeric(n_col), msg = msg)
}

.all_non_n_cols_to_char <- function(df) {
    df %>%
        dplyr::mutate(dplyr::across(!matches("^n$"), as.character))
}
