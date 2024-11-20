# Author(s): Awa Synthia
# Last modified: 2024

# Load necessary packages
library(dplyr)
library(stringr)
library(visNetwork)
library(DT)
library(plotly)

# Function to generate the InterProScan Visualization
getIPRGenesVisualization <- function(data, app_data,
                                             input_rs_iprDatabases = c("Pfam", "Phobius", "TMHMM", "Gene3D"),
                                             input_rs_iprVisType = "Analysis") {

    # Check if analysis is loaded
    if (nrow(data@df) == 0 || app_data@ipr_path == "") {
        stop("Analysis data is not loaded properly or ipr_path is missing.")
    }

    # If there is no cln_path or PcPositive column is NULL
    if (length(data@cln_path) == 0 || is.null(data@df$PcPositive)) {

        # Set column name for accessing based on the data
        n <- if ("name" %in% colnames(data)) {
            "name"
        } else {
            "AccNum"
        }
        n <- "Name"  # Hardcoded to "Name" based on original code

        # Call the `ipr2viz_web` function
        ipr_plot <- plotIPR2VizWeb(
            infile_ipr = data@ipr_path,
            accessions = data@df$Name,
            analysis = input_rs_iprDatabases,
            group_by = input_rs_iprVisType,
            name = n
        )

    } else {

        # Call the `ipr2viz` function with additional arguments
        ipr_plot <- plotIPR2Viz(
            infile_ipr = data@ipr_path,
            infile_full = data@df,
            accessions = unique(data@df$Name),
            analysis = input_rs_iprDatabases,
            group_by = input_rs_iprVisType,
            topn = 20,   # This value is hardcoded in the original code
            query = "All"
        )
    }

    # Return the plot object for further use
    return(ipr_plot)
}

# Function to generate the domain network layout visualization
getRSNetworkLayout <- function(data, app_data,
                                       cutoff = 100,
                                       layout = "nice") {

    # Check if analysis is loaded and app data has a valid ipr_path
    if (nrow(data@df) == 0 || app_data@ipr_path == "") {
        stop("Analysis not loaded or ipr_path missing.")
    }

    # Extract column names and find the first matching "DomArch" column
    cols <- colnames(data@df)
    if (is.null(cols)) {
        stop("Dataframe columns are missing.")
    }

    col <- cols[grepl("DomArch.*", cols)][1]
    if (is.null(col)) {
        stop("No domain architecture column found.")
    }

    # Clean up domain architecture columns in the data
    df_data <- data@df %>%
        mutate(across(tidyselect::starts_with("DomArch"), cleanString))

    # Generate the network using the domain_network function
    res_network <- createDomainNetwork(
        df_data,
        column = col,
        domains_of_interest = ".*",
        cutoff = cutoff,
        layout = layout
    )

    # Validate that the result is not an error
    if (any(res_network == "error")) {
      stop("Not enough nodes to construct a network.")
    }

    return(res_network)
}

# Function to generate the data table
getDataTable <- function(data) {
    if (nrow(data@df) == 0) {
        stop("No data available.
             Please ensure you have uploaded your data correctly.")
    }

    # Define the columns to be shown
    viewing_cols <- c(
        "AccNum", "QueryName", "Name", "Lineage", "Species",
        "Length", "PcPositive", "DomArch.Pfam"
    )

    # Extract the data
    d <- data@df

    # Identify columns to hide (those not in viewing_cols)
    hide_cols <- which(!colnames(d) %in% viewing_cols) - 1

    # Add hyperlinks to the "AccNum" column
    d$AccNum <- paste0("<a href='https://www.ncbi.nlm.nih.gov/protein/",
                       d$AccNum, "' target='_blank' >", d$AccNum, "</a>")

    # Create the DataTable
    dt <- DT::datatable(d,
          rownames = FALSE,
          filter = "top",  # Enable filtering
          extensions = c("Buttons"),  # Enable buttons (e.g., CSV export)
          escape = FALSE,  # Allow HTML rendering (for the hyperlinks)
          options = list(
              dom = "frlBtip",  # Layout for the table (filters, buttons, etc.)
              pageLength = 10,  # Number of rows per page
              paging = TRUE,  # Enable pagination
              # Enable regex search
              search = list(caseInsensitive = TRUE, regex = TRUE),
              language = list(
                  searchPlaceholder = "Regex or Text...",
                  filterPlaceholder = "Test"
              ),
              buttons = list(
                  list(
                      extend = "colvis",
                      text = "Add/remove column(s)"
                  ),
                  list(
                      extend = "csv",
                      text = "Download",
                      filename = "molevolvr-data",
                      exportOptions = list(
                        # Export all data, not just visible page
                          modifier = list(page = "all")
                      )
                  )
              ),
              scrollX = FALSE,  # Disable horizontal scrolling
              fixedHeader = FALSE,  # Disable fixed header
              columnDefs = list(
                # Hide columns not in 'viewing_cols'
                  list(visible = FALSE, targets = hide_cols)
              )
          )
    )

    return(dt)
}

# Function to generate query data table
getQueryDataTable <- function(query_data, query_select = NULL) {

    # Check if analysis is loaded and data is available
    if (nrow(query_data@df) == 0) {
        stop("No data available. Please ensure you have uploaded your data
             correctly.")
    }

    # Define the columns to be shown
    viewing_cols <- c("QueryName", "Species", "Lineage", "DomArch.Pfam")

    # If no specific queries are selected, display all data
    if (is.null(query_select)) {
        d <- query_data@df
    } else {
        # Filter the data for selected queries
        d <- query_data@df %>% filter(QueryName %in% query_select)
        d <- droplevels(d)
    }

    # Identify columns to hide (those not in viewing_cols)
    hide_cols <- which(!colnames(d) %in% viewing_cols) - 1

    # Add hyperlinks to the "Query" column (replace 'Query' with the actual
    # column name if different)
    d$Query <- paste0("<a href='https://www.ncbi.nlm.nih.gov/protein/",
                      d$Query, "' target='_blank' >", d$Query, "</a>")

    # Create and return the DataTable
    dt <- DT::datatable(d,
              rownames = FALSE,
              filter = "top",  # Enable filtering
              extensions = c("Buttons"),  # Enable buttons (e.g., CSV export)
              callback = DT::JS(c(
                  '$(\'div.has-feedback input[type="search"]\').attr( "placeholder", "Search..." );'
              )),
              escape = FALSE,  # Allow HTML rendering (for the hyperlinks)
              options = list(
                  dom = "frlBtip",  # Layout for the table (filters, buttons, etc.)
                  pageLength = 25,  # Number of rows per page
                  paging = TRUE,  # Enable pagination
                  # Enable regex search
                  search = list(caseInsensitive = TRUE, regex = TRUE),
                  language = list(
                      searchPlaceholder = "Regex or Text..."
                  ),
                  buttons = list(
                      list(
                          extend = "colvis",
                          text = "Add/remove column(s)"
                      ),
                      list(
                          extend = "csv",
                          text = "Download",
                          filename = "molevolvr-query",
                          exportOptions = list(
                            # Export all data, not just visible page
                              modifier = list(page = "all")
                          )
                      )
                  ),
                  scrollX = FALSE,  # Disable horizontal scrolling
                  fixedHeader = FALSE,  # Disable fixed header
                  columnDefs = list(
                    # Hide columns not in 'viewing_cols'
                      list(visible = FALSE, targets = hide_cols)
                  )
              )
    )

    return(dt)
}

# Function to read and return the FASTA file contents
readFastaData <- function(fasta_path) {

    # Check if analysis is loaded and the file path is not empty
    if (fasta_path == "" || !file.exists(fasta_path)) {
        stop("FASTA file path is invalid or the file does not exist.")
    }

    # Read the content of the FASTA file
    fasta_content <- read_file(fasta_path)

    return(fasta_content)
}

getFastaData <- function(fasta_path) {
    if (is.null(fasta_path) || fasta_path == "") {
        stop("Error: FASTA path is not provided.")
    }
    return(read_file(fasta_path))  # Read and return the FASTA data
}

# Function to get domain sequences (assumes `data` is a predefined object)
getDomData <- function(data) {
    return(data@domainSeqs)  # Return domain sequences
}

# Function to get MSA data from a given path
getMSAData <- function(msa_path) {
    if (is.null(msa_path) || msa_path == "") {
        stop("Error: MSA path is not provided.")
    }
    return(read_file(msa_path))
}

# Function to generate a heatmap
getQueryHeatmap <- function(query_data_df,
                                   heatmap_select = "All",
                                   heatmap_color = "blue") {

    # Check if analysis is loaded and query data exists
    if (nrow(query_data_df) == 0) {
        stop("No query data available.")
    }

    # Filter queries based on user selection
    if (heatmap_select == "All") {
        queries <- unique(query_data_df$QueryName)
        prot <- query_data_df %>%
            filter(Lineage != "") %>%
            tidyr::drop_na(Lineage)
    } else {
        queries <- heatmap_select
        prot <- query_data_df %>%
            filter(grepl(heatmap_select, QueryName, ignore.case = TRUE)) %>%
            filter(Lineage != "") %>%
            tidyr::drop_na(Lineage)
    }

    # Validate that the Lineage column has values
    if (all(unique(prot$Lineage) == "")) {
        stop("Lineage column not found for selected proteins.
             See the FAQ for possible reasons/solutions.")
    }

  plotLineageQuery(prot, queries = queries, colname = "QueryName",
                       cutoff = 100, color = heatmap_color)
}

# Function to retrieve domain architecture columns
getDomArchCols <- function(query_data_df) {
    # Check if query data exists
    if (nrow(query_data_df) == 0) {
        stop("No query data available.")
    }

    # Get column names
    cols <- colnames(query_data_df)

    # Filter domain architecture columns, excluding repeats
    domarch_cols <- cols[grepl("^DomArch", cols) & !grepl("repeats$", cols)]

    # Identify columns that are completely NA
    na_cols <- names(query_data_df)[apply(query_data_df, 2,
                                          function(x) all(is.na(x)))]

    # Remove NA columns from domain architecture columns
    domarch_cols <- setdiff(domarch_cols, na_cols)

    # Include SignalP_GRAM_POSITIVE if present
    if ("SignalP_GRAM_POSITIVE" %in% domarch_cols) {
        domarch_cols <- setdiff(domarch_cols, "SignalP_GRAM_POSITIVE")  # Remove first
        domarch_cols <- append(domarch_cols, "SignalP_GRAM_POSITIVE")    # Append it last
    }

    # Remove prefix from column names
    domarch_cols <- substring(domarch_cols, first = 9)

    return(domarch_cols)
}

# Function to generate main data table
getMainTable <- function(data_df, main_select = NULL) {
    # Validate input data
    if (nrow(data_df) == 0) {
        stop("No data available. Please ensure you have uploaded your data
             correctly. See Help documentation or contact JRaviLab
             (janani.ravi[AT]cuanschutz[DOT]edu).")
    }

    # Define columns to view
    viewing_cols <- c("AccNum", "QueryName", "Name", "Lineage",
                      "Species", "Length", "PcPositive",
                      "DomArch.Pfam")

    # Filter data based on selection
    if (!is.null(main_select)) {
        d <- data_df %>%
            filter(QueryName %in% main_select) %>%
            droplevels()
    } else {
        d <- data_df
    }

    # Identify columns to hide
    hide_cols <- which(!colnames(d) %in% viewing_cols)
    hide_cols <- hide_cols - 1

    # Create hyperlinks for AccNum
    d$AccNum <- paste0("<a href='https://www.ncbi.nlm.nih.gov/protein/",
                       d$AccNum, "' target='_blank'>", d$AccNum, "</a>")

    # Generate DataTable
    datatable_output <- DT::datatable(d,
                        rownames = FALSE,
                        filter = "top",
                        callback = DT::JS(c(
                            '$(\'div.has-feedback input[type="search"]\').attr( "placeholder", "Search..." );'
                        )),
                        extensions = c("Buttons"),
                        escape = FALSE,
                        options = list(
                            dom = "frlBtip",
                            pageLength = 25,
                            paging = TRUE,
                            search = list(caseInsensitive = TRUE, regex = TRUE),
                            language = list(searchPlaceholder = "Regex or Text...",
                                            filterPlaceholder = "Test"),
                            buttons = list(
                                list(extend = "colvis", text = "Add/remove column(s)"),
                                list(
                                    extend = "csv",
                                    text = "Download",
                                    filename = "molevolvr-homolog",
                                    exportOptions = list(modifier = list(page = "all"))
                                )
                            ),
                            scrollX = FALSE,
                            fixedHeader = FALSE,
                            columnDefs = list(list(visible = FALSE, targets = hide_cols))
                        )
    )

    return(datatable_output)
}

totalCounts <- function(prot, column = "DomArch", lineage_col = "Lineage",
                         cutoff = 90, RowsCutoff = FALSE, digits = 2
                         # type = "GC"
) {
  # column <- sym(column)

  prot <- select(prot, {{ column }}, {{ lineage_col }}) %>%
    filter(!is.na({{ column }}) & !is.na({{ lineage_col }})) %>%
    filter({{ column }} != "")

  prot <- summarizeByLineage(prot, column, by = lineage_col, query = "all")
  col_count <- prot %>%
    group_by(!!sym(column)) %>%
    summarise(totalcount = sum(count))

  total <- left_join(prot, col_count, by = column)

  sum_count <- sum(total$count)
  total <- total %>%
    mutate("IndividualCountPercent" = totalcount / sum_count * 100) %>%
    arrange(-totalcount, -count)

  cumm_percent <- total %>%
    select({{ column }}, totalcount) %>%
    distinct() %>%
    mutate("CumulativePercent" = 0)
  total_counter <- 0
  for (x in length(cumm_percent$totalcount):1) {
    total_counter <- total_counter + cumm_percent$totalcount[x]
    cumm_percent$CumulativePercent[x] <- total_counter / sum_count * 100
  }

  cumm_percent <- cumm_percent %>% select(CumulativePercent, {{ column }})

  total <- total %>% left_join(cumm_percent, by = as_string(column))

  # Round the percentage columns
  total$CumulativePercent <- total$CumulativePercent %>% round(digits = digits)
  total$IndividualCountPercent <- total$IndividualCountPercent %>% round(digits = digits)

  if (RowsCutoff) {
    # If total counts is being used for plotting based on number of rows,
    # don't include other observations that fall below the cummulative percent cutoff
    # , but that have the same 'totalcount' number as the cutoff observation
    total <- total %>% filter(CumulativePercent >= 100 - cutoff - .0001)
    return(total)
  }

  # Include observations that fall below the cummulative percent cutoff,
  # but that have the same 'totalcount' as the cutoff observation
  t <- total %>% filter(CumulativePercent >= 100 - cutoff)
  if (length(t) == 0) {
    cutoff_count <- 0
  } else {
    cutoff_count <- t$totalcount[nrow(t)]
  }

  total <- total %>%
    filter(totalcount >= cutoff_count) %>%
    ungroup()

  return(total)
}

getDomArchTotalCounts <- function(DA_Prot, DACutoff, DA_col, app_data) {
  # Check if ipr_path is not empty
  if (app_data@ipr_path == "") {
    stop("ipr_path is missing.")
  }

  # Calculate total counts with the specified cutoff and column
  prot_tc <- totalCounts(DA_Prot, cutoff = DACutoff, column = DA_col)

  # Replace all instances of ">" in the Lineage column with "_"
  prot_tc$Lineage <- map(prot_tc$Lineage, ~ str_replace_all(.x, ">", "_")) %>%
    unlist()

  # Return the processed data
  prot_tc
}

# Function to generate Domain Architecture Linear Table
getDomArchLinearTable <- function(DA_col, ipr_path, DA_TotalCounts_value) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Check if DA_col is provided
    if (is.null(DA_col)) {
        stop("DA_col input is required.")
    }

    da_col <- sym(DA_col)

    # Perform the data transformation
    DA_TotalCounts_value %>%
        group_by({{ da_col }}, totalcount, CumulativePercent) %>%
        summarize(LineageCount = n()) %>%
        select({{ da_col }}, LineageCount, totalcount, CumulativePercent) %>%
        arrange(-totalcount)

    # Generate the DAlin count table
    DAlin_table <- DT::datatable(
      DA_TotalCounts_value,
      selection = "single",
      extensions = c("Buttons"),
      options = list(
        pageLength = 25,
        dom = "frlBtip",
        buttons = list(
          list(
            extend = "csv",
            text = "Download",
            filename = "MolEvolvR_domarch",
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        ),
        scrollX = FALSE,
        paging = TRUE,
        fixedHeader = FALSE,
        fixedColumns = list(leftColumns = 0, rightColumns = 0)
      )
    )

    return(DAlin_table)
}

# Function to generate the Domain Architecture Lineage Plot
getDomArchHeatmapPlot <- function(DA_col, DACutoff,
                                     DA_Prot, DA_lin_color,
                                     ipr_path) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Filter the protein data for plotting
    prot <- DA_Prot %>%
        filter(Lineage != "") %>%
        drop_na(Lineage)

    # Create the plot
    plot <- plotLineageDA(
        prot,
        colname = DA_col,
        cutoff = DACutoff,
        RowsCutoff = FALSE,
        color = DA_lin_color
    )

    # Convert ggplot to plotly object
    plot <- plotly::ggplotly(plot)
    plot <- plot %>% plotly::layout(xaxis = list(side = "top"))

    return(plot)
}

# Function to generate the Domain Architecture Network
getDomNetwork <- function(DA_col, DACutoff, DA_Prot,
                                    networkLayout, ipr_path) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Prepare the selected protein data
    dn_data <- DA_Prot
    col <- sym(DA_col)  # Convert DA_col to a symbol for use with dplyr
    dn_data[[col]] <- str_replace_all(dn_data[[col]], " ", "_")

    # Filter data based on the specified column
    dn_data <- dn_data %>%
        drop_na(!!col) %>%
        filter(!!col != "")

    # Generate the domain network
    res <- createDomainNetwork(
        prot = dn_data,
        column = DA_col,
        domains_of_interest = ".*",
        cutoff = DACutoff,
        layout = networkLayout
    )

    # Validate the result
    if (any(res == "error")) {
        stop("Not enough nodes to construct a network.
             Try increasing 'Total Count Cutoff'.")
    }

    # Add export functionality
    visNetwork::addExport(res)
    res <- res %>% visExport(
        type = "png",
        name = "domain_architecture-network",
        float = "left",
        label = "Download plot",
        style = "
      background-color: #ffffff;
      color: #333333;
      border-color: #cccccc;
      padding: 6px 12px;
      font-size: 14px;
      line-height: 1.42857143;
      border-radius: 4px;
      cursor: pointer;
      display: inline-block;
      margin-bottom: 0;
      font-weight: normal;
      text-align: center;
      vertical-align: middle;
      touch-action: manipulation;
      user-select: none;
      background-image: none;
      text-decoration: none;
    "
    )

    return(res)
}

# Function to retrieve and clean Domain Architecture data
getDomArchProt <- function(app_data, DASelect) {
    # Check if the ipr_path is valid
    if (app_data@ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Validate domain architecture
    # Check if ipr_path is not empty
    if (app_data@ipr_path == "") {
        stop("ipr_path is missing.")
    }

    # Check if the data frame has rows
    if (nrow(app_data@df) <= 0) {
        stop("No data available. Please ensure you have uploaded your data correctly. See Help documentation or contact JRaviLab (janani.ravi[AT]cuanschutz[DOT]edu).")
    }

    # Check if there are domain architecture columns
    # if (length(domarch_cols) == 0) {
    #    stop("Please ensure uploaded data has domain architecture columns.")
    # }

    # Retrieve app data and clean domain architecture columns
    df_app_data <- app_data@df

    # Domain architecture column cleanup
    df_app_data <- df_app_data %>%
        mutate(across(starts_with("DomArch"), cleanString))

    # Filter based on user selection
    if (DASelect == "All") {
        return(df_app_data)
    } else {
        return(df_app_data %>% filter(grepl(DASelect, QueryName,
                                            ignore.case = TRUE)))
    }
}

# Function to retrieve domain architecture columns
getDomArchCols <- function(app_data, DASelect) {
    # Check if app data DataFrame is not empty
    if (nrow(app_data@df) <= 0) {
        stop("No data available in app data.")
    }

    # Check if ipr_path is valid
    if (app_data@ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Get column names
    cols <- colnames(app_data@df)

    # Filter domain architecture columns
    domarch_cols <- cols[grepl("^DomArch", cols) & !grepl("repeats$", cols)]

    # Retrieve the data frame
    domarch_data <- app_data@df

    # Filter data frame based on user selection
    if (DASelect != "All") {
        domarch_data <- domarch_data %>% filter(QueryName == DASelect)
    }

    # Identify columns that are completely NA or empty
    na_cols <- apply(domarch_data, 2, function(x) {
        all(is.na(x)) || all(x == "")
    })
    na_cols <- names(domarch_data)[na_cols]

    # Remove NA columns from the domain architecture columns
    domarch_cols <- setdiff(domarch_cols, na_cols)

    # Ensure "SignalP_GRAM_POSITIVE" is included if it exists
    if ("SignalP_GRAM_POSITIVE" %in% domarch_cols) {
        domarch_cols <- append(domarch_cols[!domarch_cols %in%
                                              "SignalP_GRAM_POSITIVE"],
                               "SignalP_GRAM_POSITIVE")
    }

    return(domarch_cols)
}

# Function to generate the IPR genes plot
getDomArchIPRGenesPlot <- function(app_data, da_iprDatabases,
                                       da_iprVisType, DASelect) {

    if (app_data@ipr_path == "") {
        stop("IPR path is not set.")
    }

    # Validate the data frame
    df <- app_data@df
    if (nrow(df) == 0) {
        stop("No data available. Please ensure you have uploaded
             your data correctly.")
    }

    if (is.null(da_iprDatabases) || length(da_iprDatabases) == 0) {
        stop("Please select an analysis.")
    }

    # Determine the name to use
    if (length(app_data@cln_path) == 0 || is.null(df$PcPositive)) {
        name_column <- ifelse("name" %in% colnames(df), "name", "AccNum")
        name_column <- "Name"

        # Generate the plot using the web version
        plot <- plotIPR2VizWeb(
            infile_ipr = app_data@ipr_path,
            accessions = df$Name,
            analysis = da_iprDatabases,
            group_by = da_iprVisType,
            name = name_column
        )
    } else {
        # Generate the plot using the local version
        plot <- plotIPR2Viz(
            infile_ipr = app_data@ipr_path,
            infile_full = df,
            accessions = unique(df$Name),
            analysis = da_iprDatabases,
            group_by = da_iprVisType,
            topn = 20,
            query = DASelect
        )
    }

    return(plot)
}

# Function to filter proteins for phylogeny
filterPhylogenyProteins <- function(app_data, phylo_select) {
    # Get the data frame from app_data
    df <- app_data@df

    # Check if the data frame is empty
    if (nrow(df) == 0) {
        stop("No data available in app_data.
             Please ensure it has been loaded correctly.")
    }

    # Filter the data based on user selection
    if (phylo_select == "All") {
        filtered_df <- df
    } else {
        filtered_df <- df %>%
            filter(grepl(phylo_select, QueryName, ignore.case = TRUE)) %>%
            filter(Lineage != "")
    }

    return(filtered_df)
}

# Function to retrieve representative accession numbers
getRepAccNum <- function(app_data, phylo_select,
                                                 msa_reduce_by, msa_rep_num,
                                                 rval_phylo) {

    if (rval_phylo()) {
        return(app_data@df$AccNum)
    } else {
        switch(msa_reduce_by,
               "Species" = rep_acc_species(),
               "Lineage" = rep_acc_lineage(),
               "DomArch" = {
                   if (app_data@ipr_path == "") {
                       stop("IPR path is empty. Please ensure it is set correctly.")
                   }
                   seqs <- find_top_acc(infile_full = app_data@df,
                                        n = msa_rep_num,
                                        query = phylo_select)
                   if (is.null(seqs)) {
                       stop("No sequences found.")
                   }
                   return(seqs)
               },
               "Full" = {
                   tmp <- app_data@df %>% filter(QueryName == phylo_select)
                   return(tmp$AccNum)
               },
               stop("Invalid selection for msa_reduce_by.")
        )
    }
}

getDomArchPlot <- function(ipr_path, query_names,
                                              analysis_type, group_by) {
    # Check if the input path is provided
    if (is.null(ipr_path) || ipr_path == "") {
        stop("Error: Input path for IPR data is not provided.")
    }

    # Validate that at least one domain is found
    if (length(query_domarch_cols()) < 1) {
        stop("Error: No domains found in the input sequences.")
    }

    n <- "Name"
    plot <- ipr2viz_web(
        infile_ipr = ipr_path,
        accessions = query_names,
        analysis = analysis_type,
        group_by = group_by,
        name = n
    )

    return(plot)  # Return the generated plot
}

# Function to convert accessions to names
acc2Name <- function(app_data) {
    # Check if "AccNum" is a column in the data
    if (!("AccNum" %in% colnames(app_data@df))) {
        stop("Column 'AccNum' not found in data.")
    }

    # Check if "Name" column exists and create the output data frame
    if ("Name" %in% colnames(app_data@df)) {
        # Select "AccNum" and "Name"
        df <- select(app_data@df, "AccNum", "Name")
    } else {
        # Create a data frame with "AccNum" and assign "AccNum" to "Name"
        df <- app_data@df %>%
            select("AccNum") %>%
            mutate(Name = AccNum)
    }

    return(df)
}

repAccNums <- function(phylo, msa_reduce_by, msa_rep_num, PhyloSelect, app_data) {
  # If `phylo` is true, return all `AccNum` values from `app_data`
  if (phylo) {
    return(app_data@df$AccNum)
  } else {
    # Switch based on `msa_reduce_by` value
    tmp <- filter(app_data@df)
    rep_acc_species <- createRepresentativeAccNum(tmp, reduced = "Species", accnum_col = "AccNum")

    # Get representative accession numbers by "Lineage"
    rep_acc_lineage <- createRepresentativeAccNum(tmp, reduced = "Lineage", accnum_col = "AccNum")
    switch(msa_reduce_by,
           "Species" = rep_acc_species,
           "Lineage" = rep_acc_lineage,
           "DomArch" = {
             # Check if `ipr_path` is not empty
             if (app_data@ipr_path == "") {
               stop("ipr_path is missing.")
             }
             # Find top sequences by accession number based on criteria
             seqs <- getTopAccByLinDomArch(infile_full = app_data@df, n = msa_rep_num, query = PhyloSelect)
             if (is.null(seqs)) {
               stop("No sequences found.")
             }
             return(seqs)
           },
           "Full" = {
             # Filter data for matching `QueryName` and return `AccNum`
             tmp <- app_data@df %>% filter(QueryName == PhyloSelect)
             return(tmp$AccNum)
           }
    )
  }
}

seqTree <- function(fasta_filepath){
  my_seqs <- readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
  my_seqs_msa <- msa(my_seqs)
  my_seqs_msa_aln <- msaConvert(my_seqs_msa, type="seqinr::alignment")

  #below was commented out, does it need to change as one of the parameters? the bottom keeps
  d <- dist.alignment(my_seqs_msa_aln, "identity")
  #as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

  ## Phylogenetic tree
  ## using package ape
  ## build neighbor-joining tree
  seqTree <- nj(d)
  #plot(seqTree, main="Phylogenetic Tree of MSA")
  groupInfo <- split(seqTree$tip.label,
                     gsub("_\\w+", "", seqTree$tip.label))
  seqTree <- groupOTU(seqTree, groupInfo)
  # ggtree(seqTree, aes(color=group),
  #        layout='circular') +
  #   geom_tiplab(size=1, aes(angle=angle))
  #https://yulab-smu.top/treedata-book/chapter4.html
  #offs <- 0
  tree <- ggtree(seqTree, branch.length = "dN_vs_dS") + theme_tree2(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  msaplot(tree, fasta=fasta_filepath, offset=0.5, bg_line = TRUE) + geom_tiplab(align=TRUE, linesize=0.5, size=3)

}
