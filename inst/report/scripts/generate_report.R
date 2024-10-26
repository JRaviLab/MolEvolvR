# Author(s): Awa Synthia
# Last modified: 2024

# get fasta of pathogen and/or drug
get_card_data <- function(pathogen = NULL, drug = NULL) {
  destination_dir <- "CARD_data"
  # Check if CARD data exists
  if (!dir.exists(destination_dir)) {
    dir.create(destination_dir)
  }
  if (!file.exists("CARD_data/aro_index.tsv")) {
    # Step 1: Download CARD Data
    download.file("https://card.mcmaster.ca/latest/data", "card_data.tar.bz2")
    #unzip("card_data.tar.bz2", exdir = "CARD_data")
    system(paste("tar -xjf card_data.tar.bz2 -C", destination_dir))
  }

  # Step 2: Open ARO_index.tsv
  aro_data <- read.delim("CARD_data/aro_index.tsv", sep = "\t", header = TRUE)

  # Step 3: Map CARD Short Name
  antibiotics <- read.delim("CARD_data/shortname_antibiotics.tsv",
                            sep = "\t", header = TRUE)
  pathogens <- read.delim("CARD_data/shortname_pathogens.tsv",
                          sep = "\t", header = TRUE)

  aro_data <- aro_data %>%
    mutate(
      pathogen = str_extract(CARD.Short.Name, "^[^_]+"),
      gene = str_extract(CARD.Short.Name, "(?<=_)[^_]+"),
      drug = str_extract(CARD.Short.Name, "(?<=_)[^_]+$")
    ) %>%
    left_join(pathogens, by = c("pathogen" = "Abbreviation")) %>%
    left_join(antibiotics, by = c("drug" = "AAC.Abbreviation"))

  # Sort names
  aro_data <- aro_data %>%
    arrange(Pathogen, Molecule) %>%
    group_by(Pathogen, Molecule)

  # Filter data based on user input
  filtered_data <- aro_data

  if (!is.null(pathogen)) {
    filtered_data <- filtered_data %>% filter(Pathogen == !!pathogen)
  }

  if (!is.null(drug)) {
    filtered_data <- filtered_data %>% filter(Molecule == !!drug)
  }

  # Check if filtered data is empty
  if (nrow(filtered_data) == 0) {
    stop("No data found for the specified pathogen and drug.")
  }

  # Extract protein accessions
  accessions <- filtered_data$Protein.Accession

  # Function to fetch FASTA sequence from NCBI
  get_fasta <- function(accession) {
    entrez_fetch(db = "protein", id = accession, rettype = "fasta")
  }

  # Download sequences
  fasta_sequences <- lapply(accessions, get_fasta)

  # Write sequences to a file
  writeLines(unlist(fasta_sequences), "filtered_proteins.fasta")

  return("FASTA sequences downloaded to 'filtered_proteins.fasta'.")
}

# Run analysis
run_analysis <- function(
        dupload_type = "Fasta",
        evalue = 0.00001,
        accnum_fasta_input = "",
        file_paths = list(accnum = NULL, fasta = tempfile(), msa = tempfile(),
                          blast = NULL, iprscan = NULL),
        blast_db = "refseq_protein",
        blast_hits = 100,
        blast_eval = 0.0001,
        acc_homology_analysis = TRUE,
        acc_da_analysis = TRUE,
        acc_phylogeny_analysis = FALSE,
        report_template_path = "/report/report_template.Rmd",
        output_file = file.path(tempdir(), "report.html"),
        DASelect = "All",
        mainSelect = NULL,
        PhyloSelect = NULL,
        q_heatmap_select = "All",
        DACutoff = 95,
        GCCutoff = 0.5,
        query_select = NULL,
        query_iprDatabases = NULL,
        query_iprVisType = NULL, tree_msa_tool = "ClustalO",
        levels = 2,
        DA_Col = "DomArch.Pfam",
        msa_rep_num = NULL,
        msa_reduce_by = "Species",
        rval_phylo = FALSE,
        ...

) {

    ##### Initialize Variables and Classes #####
    blast_upload_data <- new("blastUpload", df = "", seqs = "")
    ipr_upload_data <- new("iprUpload", df = "", seqs = "")
    sequence_upload_data <- new("seqUpload", seqs = "")
    data <- new("MolEvolData", msa_path = tempfile(),
                fasta_path = file_paths$fasta)
    app_data <- new("MolEvolData", msa_path = tempfile(),
                    fasta_path = file_paths$fasta)
    query_data <- new("queryData")

    if (length(dupload_type) != 1) {
        stop("dupload_type must be a single value.")
    }

    # Reset default analysis function
    resetSettings <- function() {
        acc_homology_analysis <<- FALSE
        acc_da_analysis <<- FALSE
        acc_phylogeny_analysis <<- FALSE
        domain_split <<- FALSE
    }

    fasta_set <- c("Fasta", "AccNum", "MSA")

    # Update settings based on upload type
    updateUploadType <- function(dupload_type) {
        switch(dupload_type,
               "Fasta" = {
                   resetSettings()
                   acc_homology_analysis <<- TRUE
                   acc_da_analysis <<- TRUE
               },
               "AccNum" = {
                   resetSettings()
                   acc_homology_analysis <<- TRUE
                   acc_da_analysis <<- TRUE
               },
               "MSA" = {
                   resetSettings()
                   acc_homology_analysis <<- TRUE
                   acc_da_analysis <<- TRUE
               },
               "BLAST Output" = {
                   resetSettings()
                   acc_phylogeny_analysis <<- TRUE
               },
               "InterProScan Output" = {
                   resetSettings()
                   acc_da_analysis <<- TRUE
               }
        )
    }

    updateUploadType(dupload_type)

    ####### File Upload Functions ########

    # Update function to return modified object
    uploadAccNumFile <- function(file_path, seq_obj) {
        accnum_data <- read_file(file_path)
        seq_obj@seqs <- accnum_data
        return(seq_obj)
    }

    uploadFastaFile <- function(file_path, seq_obj) {
        fasta_data <- read_file(file_path)
        seq_obj@seqs <- fasta_data
        return(seq_obj)
    }

    uploadMSAFile <- function(file_path, seq_obj) {
        msa_data <- read_file(file_path)
        seq_obj@seqs <- msa_data
        return(seq_obj)
    }

    uploadBlastFile <- function(file_path, blast_obj) {
        blast_obj@df <- file_path
        return(blast_obj)
    }

    uploadIPRScanFile <- function(file_path, ipr_obj) {
        ipr_obj@df <- file_path
        return(ipr_obj)
    }

    # Now modify the calling section to reassign the modified objects
    if (!is.null(file_paths$accnum)) {
        sequence_upload_data <- uploadAccNumFile(file_paths$accnum,
                                                 sequence_upload_data)
    }
    if (!is.null(file_paths$fasta)) {
        sequence_upload_data <- uploadFastaFile(file_paths$fasta,
                                                sequence_upload_data)
    }
    if (!is.null(file_paths$msa)) {
        sequence_upload_data <- uploadMSAFile(file_paths$msa,
                                              sequence_upload_data)
    }
    if (!is.null(file_paths$blast)) {
        blast_upload_data <- uploadBlastFile(file_paths$blast,
                                             blast_upload_data)
    }
    if (!is.null(file_paths$iprscan)) {
        ipr_upload_data <- uploadIPRScanFile(file_paths$iprscan,
                                             ipr_upload_data)
    }

    # Validation of inputs
    fasta_data <- read_file(file_paths$fasta)
    validate_and_process_inputs <- function(evalue, fasta_data) {
        if (!validate_evalue(evalue)) {
            return("Error: A numeric E-value is required. Please set a valid
                   value (e.g., 0.0001).")
        }


        if (!validate_accnum_fasta(fasta_data)) {
            return("Error: Input for AccNum/Fasta cannot be empty or invalid.")
        }

        sequence_vector <- unlist(strsplit(fasta_data, "\n"))
        return(list(message = "Inputs are valid!", sequences = sequence_vector))
    }

    validation_result <- validate_and_process_inputs(evalue, fasta_data)
    if (is.character(validation_result)) {
        return(validation_result)
    }

    # Phylogenetic Analysis Validation
    phylo <- acc_phylogeny_analysis
    if (phylo) {
        # Validate phylogenetic analysis based on upload type
        is_valid_phylo <- switch(
            dupload_type,
            "Fasta" = {
                str_count(string = sequence_upload_data@seqs, ">") > 1
            },
            "AccNum" = {
                accnums <- sequence_upload_data@seqs |>
                    strsplit("\\s*,\\s*|\\s+|\\n|,") |>
                    unlist()
                length(accnums) > 1
            },
            "MSA" = {
                str_count(string = sequence_upload_data@seqs, ">") > 1
            },
            FALSE  # Default to FALSE if no valid type is found
        )

        if (!is_valid_phylo) {
            return("Error: At least two sequences/identifiers are required for
                   phylogenetic analysis.")
        }
    }
    # Fasta-like submissions
    if (dupload_type %in% fasta_set) {
        # Can have any combination of select options
        type <- ""
        script <- ""
        postfix <- ""

        if (acc_homology_analysis && !phylo && acc_da_analysis) {
            type <- "full"
            postfix <- "full"
        }
        # Phylogenetic analysis, do full script but skip blast
        else if (phylo && acc_da_analysis) {
            type <- "phylo"
            postfix <- "phylo"
        }
        else if (acc_da_analysis) {
            type <- "da"
            postfix <- "da"
        }
        else if (acc_homology_analysis) {
            # Only run BLAST
            type <- "dblast"
            postfix <- "dblast"
        }
        else {
            # Something went wrong, throw error
          stop("Please select one of the above analyses (full, phylo, da,
               dblast) to run.")

            return()
        }
    }

    # After uploading the sequence data, you would check the uploaded data
    if (sequence_upload_data@seqs == "") {
        stop("Error: Please upload a protein sequence")
    }
    OUT_PATH <- getwd()
    unavailable_pins <- list.files(OUT_PATH)
    unavailable_pins <- strsplit_vect(unavailable_pins, pattern = "_")
    pinName <- rand_string(length = 6, post_fix = "", ignorelist = unavailable_pins)
    pin_id <- strtrim(pinName, 6)

    dir <- paste0(OUT_PATH, pin_id, "_", postfix)
    path <- paste0(dir, "/", pin_id, ".fa")
    # Adding validation logic based on upload type
    system(paste0("mkdir ", dir), wait = TRUE)
    switch(dupload_type,
           "Fasta" = {
               # Validate sequence limit
               if (str_count(sequence_upload_data@seqs, ">") > 200) {
                   stop("Error: Only submissions with less than 200 proteins
                   are accepted at this time. For analyses with more than 200 p
                        roteins please contact janani.ravi@cuanschutz.edu.")
               }

               # Validate accession numbers for FASTA submission
               if (!(validate_accnum_fasta(sequence_upload_data@seqs))) {
                   stop("Error: Please adjust the FASTA headers. Ensure a header
                        line for each sequence, no duplicate header names, and
                        no duplicate protein accession numbers.")
               }

               # TRY load fasta for fasta validation
               is_valid_aa_fasta <- tryCatch(
                   expr = {
                       tmp_file <- tempfile()
                       writeLines(sequence_upload_data@seqs, tmp_file)
                       validate_fasta(tmp_file)
                   },
                   error = function(e) {
                       warning("Error: Failed to run input FASTA verification")
                       return(FALSE)  # Return FALSE if an error occurs
                   },
                   finally = {
                       unlink(tmp_file)
                   }
               )

               # Validate AA fasta
               if (!(is_valid_aa_fasta)) {
                   stop("Error: The FASTA input could not be recognized as valid
                   Amino Acid (AA) sequence(s). MolEvolvR only accepts FASTA
                   formatted AA/protein sequences as input (not DNA/RNA).
                   For a short reference on FASTA format, review this
                   https://blast.ncbi.nlm.nih.gov/doc/blast-topics/short NCBI guide.")
               }

               write(sequence_upload_data@seqs, path)
           },
           "AccNum" = {
               # Convert Acc to FASTA
               accnum_vect <- unlist(strsplit(sequence_upload_data@seqs, "\\s*,\\s*|\\s+|\\n|,"))
               if (length(accnum_vect) > 200) {
                   stop("Error: Only submissions with less than 200 proteins are
                   accepted at this time. For analyses with more than 200
                        proteins please contact janani.ravi@cuanschutz.edu.")
               }

               # if an error is raised, validation fails
               validation_results <- purrr::map_lgl(
                   accnum_vect,
                   function(accnum) {
                       tryCatch(
                           expr = {
                               tmp <- tempfile(
                                   pattern = paste0("molevolvr_acccnum_validation-", accnum, "-", "XXXXX"),
                                   fileext = ".fa"
                               )
                               acc2fa(accnum, tmp)
                               readAAStringSet(tmp)
                               TRUE
                           },
                           error = function(e) {
                               FALSE
                           },
                           finally = {
                               suppressWarnings(sink())
                               unlink(tmp)
                           }
                       )
                   }
               )

               # If any accession numbers fail, halt submission
               if (!all(validation_results)) {
                   stop("Error: MolEvolvR could not locate sequences for the
                   following accession numbers:
                   [ {paste0(accnum_vect[which(!validation_results)], collapse = ', ')} ]
      Please try submitting FASTA sequences instead.")
               } else {
                   # Write a multifasta
                   acc2fa(accnum_vect, outpath = path)
               }
           },
           "MSA" = {
               # Validate sequence limit
               if (str_count(sequence_upload_data@seqs, ">") > 200) {
                   stop("Error: Only submissions with less than 200 proteins
                   are accepted at this time. For analyses with more than 200
                        proteins please contact janani.ravi@cuanschutz.edu.")
               }

               msa_temp <- str_replace_all(sequence_upload_data@seqs, "-", "")
               write(msa_temp, path)
           }
    )

    phylo <- if_else(phylo, "TRUE", "FALSE")

    # Clean FASTA headers
    fasta <- Biostrings::readAAStringSet(filepath = path)
    headers_original <- names(fasta)
    headers_accnum <- names(fasta) |> purrr::map_chr(function(x) extractAccNum(x))
    fasta <- cleanFAHeaders(fasta)
    headers_cleaned <- names(fasta)

    # Write a table to map original accnums to their cleaned version
    readr::write_tsv(
        x = tibble::tibble(
            "header_original" = headers_original,
            "header_accnum" = headers_accnum,
            "header_clean" = headers_cleaned
        ),
        file = file.path(dir, "query-fasta_header-map.tsv")
    )

    # Remove '*' characters which are incompatible with interproscan
    fasta <- gsub(pattern = "\\*", replacement = "X",
                  x = fasta) |> Biostrings::AAStringSet()
    Biostrings::writeXStringSet(x = fasta, filepath = path)

    if (domain_split) {
        submit_split_by_domain(
            dir = dir,
            sequences = path,
            DB = blast_db,
            NHITS = blast_hits,
            EVAL = blast_eval,
            phylo = phylo,
            type = type,
            job_code = pin_id,
        )
    } else {
        submit_full(
            dir = dir,
            sequences = path,
            DB = blast_db,
            NHITS = blast_hits,
            EVAL = blast_eval,
            phylo = phylo,
            type = type,
            job_code = pin_id,
        )
    }

    # Additional handling for different upload types
    if (dupload_type == "BLAST Output") {
        dir <- paste0(OUT_PATH, pin_id, "_blast")
        if (blast_upload_data@df == "") {
            stop("Error: Please upload a BLAST file.")
        }
        if (blast_upload_data@seqs == "" && !blast_ncbi_check) {
            stop("Error: Please provide a file containing sequences or check
                 the box to use fetch sequences for NCBI accession numbers.")
        }

        if (tools::file_ext(blast_upload_data@df) == "tsv" || blast_upload_data@df == EX_BLASTOUTPUT) {
            data <- read_tsv(blast_upload_data@df, col_names = web_blastp_hit_colnames)
        } else {
            data <- read_csv(blast_upload_data@df, col_names = web_blastp_hit_colnames)
        }

        if (nrow(data) > 5005) {
            stop("Error: BLAST output submissions are limited to 5000 proteins
            at the moment. For analyses with more than 5000 proteins please
                 contact janani.ravi@cuanschutz.edu.")
        }

        system(paste0("mkdir ", dir), wait = TRUE)
        blast_path <- paste0(dir, "/", pin_id, ".wblast.tsv")
        write_tsv(data, blast_path, col_names = FALSE)

        queries <- data$Query %>% unique()
        writeLines(queries, paste0(dir, "/", "accs.txt"))

        if (blast_ncbi_check) {
            submit_blast(
                dir = dir,
                blast = paste0(dir, "/", pin_id, ".wblast.tsv"),
                seqs = paste0(dir, "/", "seqs.fa"),
                ncbi = TRUE,
                job_code = pin_id,
                submitter_email = notify_email,
                advanced_options = isolate(rvals_advanced_options |> reactiveValuesToList() |> unlist())
            )
        } else {
            seqs <- read_file(blast_upload_data@seqs)
            writeLines(seqs, paste0(dir, "/", "seqs.fa"))
            submit_blast(
                dir = dir,
                blast = paste0(dir, "/", pin_id, ".wblast.tsv"),
                seqs = paste0(dir, "/", "seqs.fa"),
                ncbi = FALSE,
                job_code = pin_id,
                submitter_email = notify_email,
                advanced_options = isolate(rvals_advanced_options |> reactiveValuesToList() |> unlist())
            )
        }

    } else if (dupload_type == "InterProScan Output") {
        dir <- paste0(OUT_PATH, pin_id, "_ipr")
        path <- paste0(dir, "/", pin_id, "_ipr.tsv")
        if (ipr_upload_data()@df == "") {
            stop("Error: Please upload an interproscan file.")
        }
        if (ipr_upload_data()@seqs == "" && !ipr_ncbi_check) {
            stop("Error: Please provide a file containing sequences or check the
                 box to use fetch sequences for NCBI accession numbers.")
        }

        ipr <- read_tsv(ipr_upload_data()@df, col_names = FALSE)
        if (nrow(ipr) > 200) {
            stop("Error: Only submissions with less than 200 proteins are
            accepted at this time. For analyses with more than 200 proteins
                 please contact janani.ravi@cuanschutz.edu.")
        }

        system(paste0("mkdir ", dir), wait = TRUE)
        writeLines("AccNum\tSeqMD5Digest\tSLength\tAnalysis\tDB.ID\tSignDesc\tStartLoc\tStopLoc\tScore\tStatus\tRunDate\tIPRAcc\tIPRDesc\tGOTerms\tExtra\n", path)
        write_tsv(ipr, path, col_names = FALSE, append = TRUE)

        ncbi <- if_else(ipr_ncbi_check, TRUE, FALSE)
        blast <- if_else(acc_homology_analysis, TRUE, FALSE)

        if (ncbi) {
            submit_ipr(
                dir = dir,
                ipr = path,
                blast = acc_homology_analysis,
                seqs = paste0(dir, "/", "seqs.fa"),
                ncbi = TRUE,
                DB = blast_db,
                NHITS = blast_hits,
                EVAL = blast_eval,
            )
        } else {
            seqs <- read_file(ipr_upload_data()@seqs)
            writeLines(seqs, paste0(dir, "/", "seqs.fa"))
            submit_ipr(
                dir = dir,
                ncbi = FALSE,
                ipr = path,
                blast = acc_homology_analysis,
                seqs = paste0(dir, "/", "seqs.fa"),
                DB = blast_db,
                NHITS = blast_hits,
                EVAL = blast_eval,

            )
        }

    }

    # Process results for results
    fetched <- process_wrapper_dir(dir, pinName = pinName, type = type)

    # Assign fetched data to objects if the fetched list has the expected length
    if (length(fetched) == 2) {
        data <<- fetched[[1]]       # Assign the first element to 'data'
        app_data <<- fetched[[1]]   # Assign the same to 'app_data'
        query_data <<- fetched[[2]] # Assign the second element to 'query_data'

        r_nrow_initial(nrow(fetched[[1]]@df))  # Initialize row count
    }

    domarch_cols_value <- get_domarch_cols(app_data, DASelect)

    query_domarch_cols_value <- get_domarch_columns(query_data)

    mainTable_value <- generate_data_table(data)

    queryDataTable_value <- generate_query_data_table(query_data, query_select)

    fastaDataText_value <- get_fasta_data(query_data@fasta_path)

    domainDataText_value <- get_domain_data()

    msaDataText_value <- get_msa_data(query_data@msa_path)

    rs_IprGenes_value <- generate_ipr_genes_visualization(data, app_data, input_rs_iprDatabases, input_rs_iprVisType)

    rval_rs_network_layout_value <- generate_rs_network_layout(data, app_data, cutoff = 100, layout = "nice")

    rs_data_table_value <- generate_data_table(data)

    da_IprGenes_value <- generate_da_ipr_genes_plot(app_data, da_iprDatabases, da_iprVisType, DASelect)

    query_heatmap_value <- generate_query_heatmap(query_data_df, heatmap_select = "All", heatmap_color = "blue")

    DA_Prot_value <- get_DA_Prot(app_data, validate_da, DASelect)

    DALinPlot_value <- generate_DA_heatmap_plot(DA_col, DACutoff, DA_Prot, DA_lin_color, ipr_path)

    DALinTable_value <- generate_DA_lin_table(DA_col, ipr_path, DAlin_count_table_DT)

    DANetwork_value <- generate_domain_network(DA_col, DACutoff, DA_Prot, networkLayout, ipr_path)

    phylogeny_prot_value <- filter_phylogeny_proteins(app_data, phylo_select)

    acc_to_name_value <- acc_to_name(app_data)

    ####### Report Generation ########

    tryCatch({
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy(report_template_path, tempReport, overwrite = TRUE)

        # List of graphics to include in report
        params <- list(
            rs_interproscan_visualization = rs_IprGenes_value,
            proximity_network = rval_rs_network_layout_value,
            sunburst = data@df,
            data = rs_data_table_value,
            queryDataTable = queryDataTable_value,
            fastaDataText = fastaDataText_value,
            heatmap = query_heatmap_value,
            query_data = query_data,
            query_domarch_cols = query_domarch_cols_value,
            query_iprDatabases = query_iprDatabases,
            query_iprVisType = query_iprVisType,
            mainTable = maintable_value,
            DALinTable = DALinTable_value,
            DALinPlot = DALinPlot_value,
            DANetwork = DANetwork_value,
            DA_Prot = DA_Prot_value,
            domarch_cols = domarch_cols_value,
            DA_Col = DA_Col,
            DACutoff = DACutoff,
            da_interproscan_visualization = da_IprGenes_value,
            phylo_sunburst_levels = levels,
            phylo_sunburst = phylogeny_prot_value,
            tree_msa_tool = tree_msa_tool,
            rep_accnums = rep_accnums(),
            msa_rep_num = 10,
            app_data = app_data,
            PhyloSelect = PhyloSelect,
            acc_to_name = acc_to_name_value,
            rval_phylo = rval_phylo,
            msa_reduce_by = msa_reduce_by
        )

        # Render RMarkdown report
        rmarkdown::render(tempReport,
                          output_file = output_file,
                          params = params,
                          envir = new.env(parent = globalenv()))
    }, error = function(e) {
        return(paste("Error in report generation:", e$message))
    })

    return("Initialization, upload, and
           report generation completed successfully.")
}

