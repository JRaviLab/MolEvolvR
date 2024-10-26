# Author(s): Awa Synthia
# Last modified: 2024

# Load necessary packages
library(dplyr)
library(stringr)
library(visNetwork)
library(DT)
library(plotly)

# Function to generate the InterProScan Visualization
generate_ipr_genes_visualization <- function(data, app_data,
                                             input_rs_iprDatabases,
                                             input_rs_iprVisType) {

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
        ipr_plot <- ipr2viz_web(
            infile_ipr = data@ipr_path,
            accessions = data@df$Name,
            analysis = input_rs_iprDatabases,
            group_by = input_rs_iprVisType,
            name = n
        )

    } else {

        # Call the `ipr2viz` function with additional arguments
        ipr_plot <- ipr2viz(
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
generate_rs_network_layout <- function(data, app_data,
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
        mutate(across(tidyselect::starts_with("DomArch"), clean_string))

    # Generate the network using the domain_network function
    res_network <- domain_network(
        df_data,
        column = col,
        domains_of_interest = ".*",
        cutoff = cutoff,
        layout = layout
    )

    # Validate that the result is not an error
    if (res_network == "error") {
        stop("Not enough nodes to construct a network.")
    }

    return(res_network)
}

# Function to generate the data table
generate_data_table <- function(data) {
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
generate_query_data_table <- function(query_data, query_select = NULL) {

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
read_fasta_data <- function(fasta_path) {

    # Check if analysis is loaded and the file path is not empty
    if (fasta_path == "" || !file.exists(fasta_path)) {
        stop("FASTA file path is invalid or the file does not exist.")
    }

    # Read the content of the FASTA file
    fasta_content <- read_file(fasta_path)

    return(fasta_content)
}

get_fasta_data <- function(fasta_path) {
    if (is.null(fasta_path) || fasta_path == "") {
        stop("Error: FASTA path is not provided.")
    }
    return(read_file(fasta_path))  # Read and return the FASTA data
}

# Function to get domain sequences (assumes `data` is a predefined object)
get_domain_data <- function() {
    return(data@domainSeqs)  # Return domain sequences
}

# Function to get MSA data from a given path
get_msa_data <- function(msa_path) {
    if (is.null(msa_path) || msa_path == "") {
        stop("Error: MSA path is not provided.")
    }
    return(read_file(msa_path))
}

# Function to generate a heatmap
generate_query_heatmap <- function(query_data_df,
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

    # Assuming `lineage.Query.plot` is a custom function for plotting
    lineage.Query.plot(prot, queries = queries, colname = "QueryName",
                       cutoff = 100, color = heatmap_color)
}

# Function to retrieve domain architecture columns
get_domarch_columns <- function(query_data_df) {
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
generate_main_table <- function(data_df, main_select = NULL) {
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

# Function to generate Domain Architecture Linear Table
generate_DA_lin_table <- function(DA_col, ipr_path, DAlin_count_table_DT) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Check if DA_col is provided
    if (is.null(DA_col)) {
        stop("DA_col input is required.")
    }

    # Generate the DAlin count table
    DAlin_table <- DAlin_count_table_DT()

    return(DAlin_table)
}

# Function to generate the Domain Architecture Lineage Plot
generate_DA_heatmap_plot <- function(DA_col, DACutoff,
                                     DA_Prot, DA_lin_color,
                                     analysis_loaded, ipr_path) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Filter the protein data for plotting
    prot <- DA_Prot() %>%
        filter(Lineage != "") %>%
        drop_na(Lineage)

    # Create the plot
    plot <- lineage.DA.plot(
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
generate_domain_network <- function(DA_col, DACutoff, DA_Prot,
                                    networkLayout, ipr_path) {
    # Check if ipr_path is valid
    if (ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Prepare the selected protein data
    dn_data <- DA_Prot()
    col <- sym(DA_col)  # Convert DA_col to a symbol for use with dplyr
    dn_data[[col]] <- str_replace_all(dn_data[[col]], " ", "_")

    # Filter data based on the specified column
    dn_data <- dn_data %>%
        drop_na(!!col) %>%
        filter(!!col != "")

    # Generate the domain network
    res <- domain_network(
        prot = dn_data,
        column = col,
        domains_of_interest = ".*",
        cutoff = DACutoff,
        layout = networkLayout
    )

    # Validate the result
    if (res == "error") {
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
get_DA_Prot <- function(app_data, validate_da, DASelect) {
    # Check if the ipr_path is valid
    if (app_data@ipr_path == "") {
        stop("InterPro path is empty.")
    }

    # Validate domain architecture
    validate_da()

    # Retrieve app data and clean domain architecture columns
    df_app_data <- app_data@df

    # Domain architecture column cleanup
    df_app_data <- df_app_data %>%
        mutate(across(starts_with("DomArch"), clean_string))

    # Filter based on user selection
    if (DASelect == "All") {
        return(df_app_data)
    } else {
        return(df_app_data %>% filter(grepl(DASelect, QueryName,
                                            ignore.case = TRUE)))
    }
}

# Function to retrieve domain architecture columns
get_domarch_cols <- function(app_data, DASelect) {
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
generate_da_ipr_genes_plot <- function(app_data, da_iprDatabases,
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
        plot <- ipr2viz_web(
            infile_ipr = app_data@ipr_path,
            accessions = df$Name,
            analysis = da_iprDatabases,
            group_by = da_iprVisType,
            name = name_column
        )
    } else {
        # Generate the plot using the local version
        plot <- ipr2viz(
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
filter_phylogeny_proteins <- function(app_data, phylo_select) {
    # Validate the input app_data
    if (!analysis_loaded()) {
        stop("Analysis not loaded.")
    }

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
get_representative_accession_numbers <- function(app_data, phylo_select,
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

generate_domain_architecture_plot <- function(ipr_path, query_names,
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
acc_to_name <- function(app_data) {
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

