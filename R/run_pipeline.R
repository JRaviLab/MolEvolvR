# Author(s): Awa Synthia
# Last modified: 2024

# Load necessary libraries
library(httr)
library(data.table)
library(readr)
library(rentrez)

getSeqs <- function(sequences,
                          acc_file_path = "accs.txt",
                          dir_path = "~",
                          separate = TRUE) {
  seqs <- readAAStringSet(sequences)
  cln_names <- c()
  for (accnum in names(seqs)) {
    if (grepl("\\|", accnum)) {
      accnum_cln <- strsplit(accnum, "\\|")[[1]][2]
      accnum_cln <- strsplit(accnum_cln, " ")[[1]]
    } else {
      accnum_cln <- strsplit(accnum, " ")[[1]][1]
    }
    cln_names <- append(cln_names, accnum_cln)
    write(accnum_cln, file = acc_file_path, append = TRUE)
    if (separate) {
      write(paste0(dir_path, "/", accnum_cln, ".faa"),
            file = "input.txt", append = TRUE)
      write(paste0(">", accnum_cln),
            file = paste0(accnum_cln, ".faa"), append = TRUE)
      write(toString(seqs[accnum]),
            file = paste0(accnum_cln, ".faa"), append = TRUE)
    }
  }
  names(seqs) <- cln_names
  writeXStringSet(seqs, sequences, format = "fasta")
  return(length(seqs))
}
runFull <- function(
    dir = "/data/scratch",
    DB = Sys.getenv("BLAST_DB", unset = "refseq"),
    NHITS = Sys.getenv("BLAST_HITS", unset = 100),
    EVAL = Sys.getenv("BLAST_EVALUE", unset = 0.00001),
    sequences = "~/test.fa",
    phylo = "FALSE",
    by_domain = "FALSE",
    domain_starting = "~/domain_seqs.fa",
    type = "full",
    job_code=NULL,
    submitter_email=NULL,
    advanced_options=NULL,
    get_slurm_mails=FALSE
) {
  # Set working directory
  setwd(dir)

  advanced_options_names <- names(advanced_options[advanced_options == TRUE])

  # Write job submission params to file
  job_args <- list(
    submission_type = type,
    database = ifelse(phylo == FALSE, DB, NA),
    nhits = ifelse(phylo == FALSE, NHITS, NA),
    evalue = ifelse(phylo == FALSE, EVAL, NA),
    submitter_email = submitter_email,
    advanced_options = advanced_options_names,
    job_code = job_code
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  # Create a log file
  write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa
        \tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2DomArch\tduration",
        "logfile.tsv")

  # Process sequences (local handling)
  if (phylo == "FALSE") {
    # Split the sequences if needed, store them locally
    num_seqs <- getSeqs(sequences, dir_path = dir, separate = TRUE)

    fasta <- Biostrings::readAAStringSet(sequences)
    headers_original <- names(fasta)
    headers_accnum <- names(fasta) |> purrr::map_chr(function(x) extractAccNum(x))

    # Execute BLAST locally (instead of submitting jobs to the cluster)
    for (i in 1:num_seqs) {
      # Assume each sequence is saved separately
      input_file <- paste0(dir, "/", headers_accnum[i], ".faa")
      # output_file <- paste0(dir, "/blast_output_", i, ".txt")

      # Construct the local BLAST command (make sure 'blastn' is available locally)
      runMolevolvrPipeline(input_file, DB, NHITS, EVAL, is_query = F, type, i)

      #cmd <- sprintf(
      #  "deltablast -query %s -db %s -out %s -num_alignments %d -evalue %f -remote",
      #  input_file, DB, output_file, NHITS, EVAL
      #)

      # Execute BLAST locally
      #system(cmd)

      cat(sprintf("BLAST for sequence %d completed.\n", i), file=stderr())
    }
  } else {
    # Handle phylogenetic analysis if needed
    cat("Phylogenetic analysis is not supported in the local version yet.\n")
  }

  # Simulate query run locally
  runMolevolvrPipeline(sequences, DB, NHITS, EVAL, is_query = TRUE, type)
  # cmd_query <- sprintf(
  #   "deltablast -query %s -db %s -out %s_query.txt -num_alignments %d -evalue %f -remote",
  #   sequences, DB, paste0(dir, "/query_output"), NHITS, EVAL
  # )
  # system(cmd_query)

  cat("Query analysis completed.\n")

  # Status update
  num_runs <- num_seqs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")

  cat("All analyses completed locally.\n")
}

# Define the main pipeline function
runMolevolvrPipeline <- function(input_paths, db, nhits, eval,
                                   is_query, type, i) {

  # Start time
  start <- Sys.time()
  OUTPATH <- getwd()   # Set output directory to the current working directory

  # If IS_QUERY is True, handle query data
  if (is_query == TRUE) {

    FILE <- input_paths
    PREFIX <- "query_data"
    OUTDIR <- file.path(OUTPATH, PREFIX)
    dir.create(OUTDIR, showWarnings = FALSE)  # Create output directory

    all_accnums_file <- file.path(OUTDIR, paste0(PREFIX, ".all_accnums.fa"))
    file.copy(FILE, all_accnums_file)

    # Parse accession numbers
    accnums_file <- file.path(OUTPATH, "query-fasta_header-map.tsv")
    parsed_accnums_file <- file.path(OUTDIR, "parsed_accnums.txt")

    if (file.exists(accnums_file)) {
      parsed_accnums <- read.table(accnums_file,
                                   header = FALSE, sep = "\t",
                                   stringsAsFactors = FALSE)
      writeLines(parsed_accnums$V2[-1], parsed_accnums_file)  # Skip header
    }

    # Copy starting_accs.txt if exists, or fallback to accs.txt
    if (file.exists("starting_accs.txt")) {
      file.copy("starting_accs.txt",
                file.path(OUTDIR, paste0(PREFIX, ".all_accnums.txt")))
    } else {
      file.copy("accs.txt",
                file.path(OUTDIR, paste0(PREFIX, ".all_accnums.txt")))
      file.copy(parsed_accnums_file,
                file.path(OUTDIR, paste0(PREFIX, ".parsed_accnums.txt")))
    }

    # setwd(OUTDIR)

    # Run acc2info
    runAcc2Info(parsed_accnums_file, PREFIX, OUTDIR)

    replaceAccNums(file.path(OUTDIR, paste0(PREFIX, ".acc2info.tsv")),
                              file.path(OUTPATH, "query-fasta_header-map.tsv"),
                              file.path(OUTDIR, paste0(PREFIX, ".acc2info.tsv")))

    file.copy(file.path(OUTDIR, paste0(PREFIX, ".acc2info.tsv")),
              file.path(OUTDIR, paste0(PREFIX, ".blast.cln.tsv")))

  } else {
    # Handle homolog data
    FILE <- readLines(input_paths)[1]  # Assume single file for local usage
    F_value <- basename(FILE)
    PREFIX <- sub("\\.faa$", "", basename(FILE))
    PREFIX <- gsub(">", "", PREFIX)  # Extract prefix from file name
    OUTDIR <- file.path(OUTPATH, paste0(PREFIX))
    dir.create(OUTDIR, showWarnings = FALSE)  # Create output directory
    # setwd(OUTDIR)

    # Run DELTABLAST
    runDELTABLAST(input_paths, PREFIX, OUTDIR, db, nhits, eval)

    # Run ACC2FA
    convertAccNum2Fasta(file.path(OUTDIR, paste0(PREFIX, ".dblast.tsv")),
                            PREFIX, OUTDIR)

    # Run ACC2INFO
    runAcc2Info(file.path(OUTDIR, paste0(PREFIX, ".all_accnums.txt")),
                 PREFIX, OUTDIR)

    # Clean up BLAST results
    cleanupBlast(file.path(OUTDIR, paste0(PREFIX, ".dblast.tsv")),
                  file.path(OUTDIR, paste0(PREFIX, ".acc2info.tsv")),
                  PREFIX, F)

  }

  # Sys.sleep(30)

  # Run BLASTCLUST
  runCDHIT(file.path(OUTDIR, paste0(PREFIX, ".all_accnums.fa")),
                 PREFIX, OUTDIR )

  # Convert clusters to table
  clust2Table(file.path(OUTDIR, paste0(PREFIX, ".bclust.L60S80.tsv")),
            file.path(OUTDIR, paste0(PREFIX, ".blast.cln.tsv")))

  # Run INTERPROSCAN
  runIPRScan2(file.path(OUTDIR, paste0(PREFIX, ".all_accnums.fa")),
                   PREFIX, OUTDIR)
  new_header <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis", "DB.ID",
                  "SignDesc", "StartLoc", "StopLoc", "Score",
                  "Status", "RunDate", "IPRAcc", "IPRDesc")

  temp_data <- read_tsv(file.path(OUTDIR, paste0(PREFIX, ".iprscan.tsv")),
                        col_names = FALSE)

  colnames(temp_data) <- new_header

  write.table(temp_data, file.path(OUTDIR, paste0(PREFIX, ".iprscan.tsv")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  # Run IPR2LIN
  ipr2Linear(file.path(OUTDIR, paste0(PREFIX, ".iprscan.tsv")),
          file.path(OUTDIR, paste0(PREFIX, ".acc2info.tsv")), PREFIX)


  # Optionally run ipr2DomArch if IS_QUERY is true
  if (is_query == TRUE) {
    ## perform ipr2DomArch on iprscan results
    da <- ipr2DomArch(file.path(OUTDIR, paste0(PREFIX, ".iprscan_cln.tsv")),
                 PREFIX, "NA")

    ## if blast results are provided, call appendIPR
    if (is.null(file.path(OUTDIR, paste0(PREFIX, ".iprscan_cln.tsv"))) |
        is.na(PREFIX)) {
      print("No blast results provided, moving on.")
    } else {
      file.copy(file.path(OUTDIR, paste0(PREFIX, ".ipr_domarch.tsv")),
                file.path(OUTDIR, paste0(PREFIX, ".full_analysis.tsv")))
    }
  } else {
    # perform ipr2DomArch on iprscan results
    da <- ipr2DomArch(file.path(OUTDIR, paste0(PREFIX, ".iprscan_cln.tsv")),
                 PREFIX)

    ## if blast results are provided, call appendIPR
    if (is.null(file.path(OUTDIR, paste0(PREFIX, ".iprscan_cln.tsv"))) |
        is.na(PREFIX)) {
      print("No blast results provided, moving on.")
    } else {
      appendIPR(ipr_da = da,
                 blast =  file.path(OUTDIR, paste0(PREFIX, ".cln.clust.tsv")),
                 prefix = PREFIX)
    }
  }

  # Copy the input fasta file to the output directory
  # file.copy(input_paths, OUTDIR)

  # Total run time
  dur <- difftime(Sys.time(), start, units = "secs")
  cat("\nTotal run time:", dur, "seconds\n")

  # Log the run times
  logfile_path <- file.path(OUTPATH, "logfile.tsv")
  log_data <- data.frame(
    START_DT = format(start, "%d/%m/%Y-%H:%M:%S"),
    STOP_DT = format(Sys.time(), "%d/%m/%Y-%H:%M:%S"),
    PREFIX = PREFIX,
    dur = dur
  )

  write.table(log_data, file = logfile_path, sep = "\t",
              row.names = FALSE, col.names = FALSE,
              append = TRUE, quote = FALSE)
}

# Define the acc2info function
acc2info <- function(infile, prefix, outdir) {
  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  outfile <- file.path(outdir, paste0(prefix, ".acc2info.tsv"))

  # Print column names to the output file
  # Print column names to outfile

  # Read the input file of accession numbers
  acc_nums <- readLines(infile)

  # Random sleep to avoid overloading the server
  # Sys.sleep(sample(1:10, 1))

  # EPOST to get WebEnv and QueryKey
  epost_response <- entrez_post(db = "protein", id = acc_nums)
  webenv <- epost_response$WebEnv
  query_key <- epost_response$QueryKey
  # EFETCH to get the document summaries
  docsums <- entrez_summary(db = "protein", web_history = epost_response)

  # Check if any atomic values exist in docsums
  any_atomic <- FALSE
  # Check if any value in the current docsum is atomic
  for (docsum in docsums) {
    # Check if any value in the current docsum is atomic
    if (is.atomic(docsum)) {
      any_atomic <- TRUE
      break  # Exit loop early if an atomic value is found
    }
  }


  if (any_atomic) {
    parsed_data <- data.frame(
      AccNum = docsums$oslt$value,
      AccNum.noV = docsums$caption,
      FullAccNum = docsums$extra,
      Description = docsums$title,
      Length = docsums$slen,
      TaxID = docsums$taxid,
      Species = docsums$organism,
      SourceDB = ifelse(is.null(docsums$sourcedb), NA,
                        docsums$sourcedb),
      Completeness = ifelse(is.null(docsums$completeness), NA,
                            docsums$completeness),
      stringsAsFactors = FALSE  # Avoid factors in data frame
    )
  } else {
    # If no atomic values exist, use lapply to parse each docsum
    parsed_data <- do.call(rbind, lapply(docsums, function(docsum) {
      tryCatch({
        data.frame(
          AccNum = docsum$oslt$value,
          AccNum.noV = docsum$caption,
          FullAccNum = docsum$extra,
          Description = docsum$title,
          Length = docsum$slen,
          TaxID = docsum$taxid,
          Species = docsum$organism,
          SourceDB = ifelse(is.null(docsum$sourcedb), NA,
                            docsum$sourcedb),
          Completeness = ifelse(is.null(docsum$completeness), NA,
                                docsum$completeness),
          stringsAsFactors = FALSE  # Avoid factors in data frame
        )
      }, error = function(e) {
        # Return a data frame filled with NA values in case of an error
        return(data.frame(
          AccNum = NA,
          AccNum.noV = NA,
          FullAccNum = NA,
          Description = NA,
          Length = NA,
          TaxID = NA,
          Species = NA,
          SourceDB = NA,
          Completeness = NA,
          stringsAsFactors = FALSE
        ))
      })
    }))
  }

  # Parse the fetched data
  # Check if we got any data
  if (nrow(parsed_data) == 0) {
    cat("No data found for accession numbers in NCBI. Trying UniProt...\n")

    # Split input file into smaller chunks
    split_files <- split(acc_nums, ceiling(seq_along(acc_nums) / 100))

    for (i in seq_along(split_files)) {
      accnum <- paste(split_files[[i]], collapse = ",")  # Join with comma
      url <- paste0("https://www.ebi.ac.uk/proteins/api/proteins?accession=",
                    accnum)

      response <- GET(url, accept("application/xml"))

      if (status_code(response) == 200) {
        xml_content <- content(response, as = "text")

        # Parse the XML response
        parsed_uniprot <- xml2::read_xml(xml_content)

        # Extract required elements (customize based on XML structure)
        entries <- xml2::xml_find_all(parsed_uniprot, "//entry")

        for (entry in entries) {
          accession <- xml2::xml_text(xml2::xml_find_first(entry,
                                                           ".//accession"))
          full_name <- xml2::xml_text(xml2::xml_find_first(entry,
                                                           ".//fullName"))
          length_seq <- xml2::xml_attr(xml2::xml_find_first(entry,
                                                            ".//sequence"),
                                                            "length")
          db_reference <- xml2::xml_attr(xml2::xml_find_first(entry,
                                                              ".//dbReference"),
                                                              "id")
          dataset <- xml2::xml_attr(entry, "dataset")
          name <- xml2::xml_text(xml2::xml_find_first(entry, ".//name"))

          # Create a row of data
          new_row <- data.frame(AccNum = accession,
                                AccNum.noV = gsub("\\|", "", accession),
                                FullAccNum = accession,
                                Description = full_name,
                                Length = as.integer(length_seq),
                                TaxID = db_reference,
                                Species = name,
                                SourceDB = dataset,
                                Completeness = "NA",  # Adjust as needed
                                stringsAsFactors = FALSE)

          # Append to the outfile
          write.table(new_row, file = outfile, sep = "\t", col.names = TRUE,
                      row.names = FALSE, quote = FALSE, append = TRUE)
        }
      }
    }
  } else {
    # Append NCBI data to outfile
    write.table(parsed_data, file = outfile, sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE, append = TRUE)
  }

  cat("Data saved to:", outfile, "\n")
}

acc2InfoPhylo <- function(infile, outdir) {
  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  outfile <- file.path(outdir, "acc2info.tsv")

  # Create a data frame to store results
  results <- data.frame(
    Caption = character(),
    Extra = character(),
    Title = character(),
    Length = character(),
    TaxID = character(),
    Organism = character(),
    SourceDB = character(),
    Completeness = character(),
    stringsAsFactors = FALSE
  )

  # Process accession numbers in batches to avoid overloading
  for (acc in acc_nums) {
    # Fetch data using epost and efetch
    response <- GET(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&format=xml&id=%s",
                            acc))

    # Check if the request was successful
    if (response$status_code == 200) {
      # Parse the XML response
      parsed_xml <- read_xml(content(response, "text"))

      # Extract relevant information
      doc_summary <- xml_find_all(parsed_xml, ".//DocumentSummary")

      for (summary in doc_summary) {
        results <- rbind(results, data.frame(
          Caption = xml_text(xml_find_first(summary, ".//Caption")),
          Extra = xml_text(xml_find_first(summary, ".//Extra")),
          Title = xml_text(xml_find_first(summary, ".//Title")),
          Length = xml_text(xml_find_first(summary, ".//Slen")),
          TaxID = xml_text(xml_find_first(summary, ".//TaxId")),
          Organism = xml_text(xml_find_first(summary, ".//Organism")),
          SourceDB = xml_text(xml_find_first(summary, ".//SourceDb")),
          Completeness = xml_text(xml_find_first(summary, ".//Completeness")),
          stringsAsFactors = FALSE
        ))
      }
    } else {
      warning(sprintf("Failed to fetch data for accession %s: %s",
                      acc, response$status_code))
    }
  }

  # Write results to the output file
  write.table(results, file = outfile, sep = "\t", row.names = FALSE,
              quote = FALSE)
}

# Main function to run based on the prefix
runAcc2Info <- function(infile, prefix, outdir) {
  if (prefix == "NA") {
    acc2InfoPhylo(infile, outdir)
  } else {
    acc2info(infile, prefix, outdir)
  }
}

subsAccnum4cc2Info <- function(df_acc2info, df_header_map) {
  df_result <- df_header_map |>
    # set column name in header map to match accnum col in acc2info
    dplyr::rename(AccNum = header_accnum) |>
    # join onto header map
    dplyr::left_join(df_acc2info, by = "AccNum") |>
    # deselect accnum
    dplyr::select(-AccNum) |>
    # set the accnum col to the cleaned form
    dplyr::rename(AccNum = header_clean) |>
    # rm excess columns from header map file
    dplyr::select(-header_original)
  return(df_result)
}


replaceAccNums <- function(path_acc2info,
                                      path_query_header_map, path_out) {

  # Read the input files
  df_acc2info <- read_tsv(path_acc2info)
  df_query_header_map <- read_tsv(path_query_header_map)

  # Substitute accession numbers
  df_acc2info_substituted <- subsAccnum4cc2Info(df_acc2info,
                                                            df_query_header_map)

  # Print the substituted dataframe
  print("### df_acc2info_substituted")
  print(df_acc2info_substituted)

  # Write the substituted dataframe to the output file
  write_tsv(df_acc2info_substituted, file = path_out, col_names = TRUE)
}


runDELTABLAST <- function(infile, prefix, outdir,
                           db = "refseq_protein",
                           nhits = 5000, evalue = 1e-5,
                           threads = 10) {

  # Prepare output file path
  # outfile <- file.path(outdir, paste0(prefix, ".dblast.tsv"))

  outfile <- paste0(outdir, "/" , prefix, ".dblast.tsv")

  # Print I/O messages
  cat("\nNow processing:", infile, "\n")
  cat("Running against:", db, "\n")
  cat("E-value: ≤", evalue, "\n")
  cat("Top", nhits, "hits/alignments\n")
  cat("Output filepath:", outfile, "\n")

  # Designating database based on user input
  dblast_db <- switch(db,
                      nr = "nr",
                      refseq_protein = "refseq_protein",
                      stop("Invalid database specified.
                           Choose either 'nr' or 'refseq'."))

  # Core command for running DELTABLAST
  command <- sprintf(
    "deltablast -query %s -db %s -out %s -num_alignments %d -evalue %s -outfmt '6 qacc sacc sseqid sallseqid stitle sscinames staxids pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore ppos' -remote",
    infile, dblast_db, outfile, nhits, evalue
  )
  # Run the DELTABLAST command
  system(command, intern = TRUE)

  cat("DELTABLAST completed.\n")
}


# This script converts AccNum to Fasta using NCBI's EDirect or EBI's API

convertAccNum2Fasta <- function(infile, prefix, outdir) {

  # Create the output file path
  outfile <- file.path(outdir, paste0(prefix, ".all_accnums.fa"))

  # Print status
  cat("\n####################\n")
  cat("BEGIN EDIRECT SEARCH\n")
  cat("####################\n")
  cat("Processing input file:", infile, "\n")

  # Create temporary files based on the number of columns in input file
  temp_acc_file <- file.path(outdir, paste0(prefix, ".all_accnums.txt"))
  cols <- length(strsplit(readLines(infile, n = 1), "\t")[[1]])

  if (cols > 1) {
    # Extract unique homolog accession numbers from the 2nd column
    cat("Creating temp file with 2nd column unique accessions...\n")
    system(sprintf("awk -F '\\t' '{ print $2 }' %s | sort -u > %s", infile, temp_acc_file))
  } else {
    # Extract unique homolog accession numbers from the 1st column
    cat("Creating temp file with 1st column unique accessions...\n")
    system(sprintf("awk -F '\\t' '{ print $1 }' %s | sort -u > %s", infile, temp_acc_file))
  }

  # Split accessions into chunks of 1000
  system(sprintf("split -l 1000 -e %s %s/acc", temp_acc_file, outdir))

  # Fetch FASTA sequences
  cat("\nObtaining FASTA files\n")
  acc_files <- list.files(path = outdir, pattern = "^acc", full.names = TRUE)

  for (x in acc_files) {
    # Accession numbers in a comma-separated format
    accnum <- paste(readLines(x), collapse = ",")
    fasta_sequence <- entrez_fetch(db = "protein", id = accnum,
                                   rettype = "fasta")
    # Write the fetched sequence to the output file
    writeLines(fasta_sequence, con = outfile)
  }

  # Check if any sequences were retrieved
  num_seqs <- as.numeric(system(sprintf("grep '>' %s | wc -l", outfile),
                                intern = TRUE))

  # If no sequences found, try fetching from EBI API
  if (num_seqs < 1) {
    cat("\nNo sequences retrieved. Trying EBI API...\n")
    for (x in acc_files) {
      accnum <- paste(readLines(x), collapse = ",")
      system(sprintf(
        "curl -X GET --header 'Accept:text/x-fasta' 'https://www.ebi.ac.uk/proteins/api/proteins?accession=%s' >> %s",
        accnum, outfile
      ))

    }
  }

  # Remove temporary files
  cat("\nRemoving temporary files\n")
  system(sprintf("rm %s/acc*", outdir))

  # Print end message
  cat("#####################\n")
  cat("END OF EDIRECT SEARCH\n")
  cat("#####################\n")
}


cleanupBlast <- function(infile_blast, acc2info, prefix, wblast = F) {

  outdir <- dirname(infile_blast)
  # Load and clean acc2info file
  acc2info_out <- fread(input = acc2info, sep = "\t",
                        header = TRUE, fill = TRUE) %>%
    mutate(FullAccNum = gsub("\\|", "", FullAccNum)) %>%
    mutate(FullAccNum = gsub(".*[a-z]", "", FullAccNum))

  # If using web-blast output
  if (wblast == TRUE) {
    blast_out <- fread(input = infile_blast, sep = "\t", header = FALSE,
                       col.names = web_blastp_hit_colnames, fill = TRUE)
    cleanedup_blast <- blast_out %>%
      mutate(AccNum = gsub("\\|", "", AccNum)) %>%
      mutate(AccNum = gsub(".*[a-z]", "", AccNum)) %>%
      mutate(PcIdentity = round(as.double(PcIdentity), 2))

    # Merge cleaned blast output with acc2info
    cleanedup_blast <- merge(cleanedup_blast, acc2info_out,
                             by.x = "AccNum",
                             by.y = "FullAccNum",
                             all.x = TRUE)
    names(cleanedup_blast)[names(cleanedup_blast) == "Species.y"] <- "Species"

    # Additional calculations for PcPositive
    cleanedup_blast <- cleanedup_blast %>%
      mutate(PcPosOrig = as.numeric(PcPosOrig)) %>%
      mutate(AlnLength = as.numeric(AlnLength)) %>%
      mutate(PcPositive = PcPosOrig)

  } else if (wblast == FALSE) {
    # If using classic-blast output
    blast_out <- read_tsv(file = infile_blast, col_names = cl_blast_colnames)
    cleanedup_blast <- blast_out %>%
      mutate(AccNum = gsub("\\|", "", AccNum)) %>%
      mutate(AccNum = gsub(".*[a-z]", "", AccNum)) %>%
      mutate(Species = gsub(";.*$", "", Species)) %>%
      mutate(PcIdentity = round(PcIdentity, 2)) %>%
      mutate(PcPositive = round((PcPosOrig * AlnLength / QLength), digits = 2))

    # Merge cleaned blast output with acc2info
    cleanedup_blast <- merge(cleanedup_blast, acc2info_out,
                             by.x = "AccNum",
                             by.y = "FullAccNum",
                             all.x = TRUE) %>%
      select(-Species.x, -TaxID.x)
    names(cleanedup_blast)[names(cleanedup_blast) == "Species.y"] <- "Species"
  }

  # TaxID to lineage mapping
  cleanedup_blast$TaxID <- as.integer(cleanedup_blast$TaxID)
  lineage_map <- fread(
    system.file("common_data", "lineage_lookup.txt", package = "MolEvolvR", mustWork = TRUE),
    header = TRUE, fill = TRUE, colClasses = lineage_map_cols)

  # Merge with lineage map and clean up columns
  mergedLins <- merge(cleanedup_blast, lineage_map, by = "TaxID",
                      all.x = TRUE) %>%
    mutate(Species = Species.y, Spp.blast = Species.x) %>%
    select(any_of(cl_blast_postcln_cols))

  # Add names and prepare the output
  blast_names <- addName(mergedLins)

  # Create output file name
  file_name <- file.path(outdir, paste0(prefix, ".blast.cln.tsv"))

  # Write cleaned data to file
  write_tsv(blast_names, file_name, col_names = TRUE)
}


# Function to run BLASTCLUST on given input
runCDHIT <- function(infile, suffix, outdir) {

  # Prepare output file path
  outfile <- file.path(outdir, paste0(suffix, ".bclust.L60S80.tsv"))
  input_file <- file.path(outdir, paste0(suffix, ".bclust.L60S80.tsv.clstr"))


  # Print process messages
  cat("\n#####################################\n")
  cat("## Now running BLASTCLUST on file(s):", infile, "\n")
  cat("#####################################\n")

  # Cluster sequences w/ BLASTCLUST/CD-HIT
  # blastclust_cmd <- paste("blastclust -i", infile, "-o", outfile, "-p T -L .6 -b T -S 80 -a 8")
  cdhit_command <- sprintf(
    "cd-hit -i %s -o %s -c 0.8 -aS 0.6 -T 8",
    infile, outfile
  )
  # clean_cdhit_format <- sprintf(
  #   "awk '/^>Cluster/ {if(NR>1)printf \"\\n\"; next} /WP_/ {start=index($0, \"WP_\"); if(start) {end=index(substr($0, start), \"...\"); if (end == 0) end=length($0); printf \"%%s \", substr($0, start, end-1)}}' %s > %s",
  #   input_file,
  #   outfile
  # )

  cat("\nPerforming BLASTCLUST analysis on", infile, "\n")
  # system(blastclust_cmd)

  system(cdhit_command)
  # system(clean_cdhit_format)
  cleanCDHIT(input_file, outfile)
  cat("Cleaning CDHIT results completed")
}

# extract_sequences("input_file.tsv.clstr", "output_file.tsv")
cleanCDHIT <- function(input_file, output_file) {
  # Read the input file line by line
  lines <- readLines(input_file)

  # Initialize an empty vector to store extracted identifiers
  extracted <- c()

  # Loop through lines and extract content between ">" and "..."
  for (line in lines) {
    # Only process lines containing valid sequences (skip cluster headers)
    if (grepl(">", line) && grepl("\\.\\.\\.", line)) {
      # Extract the sequence identifier between ">" and "..."
      match <- regmatches(line, regexpr(">(.*?)\\.\\.\\.", line))
      if (length(match) > 0) {
        sequence <- gsub("^>", "", match)  # Remove the leading ">"
        sequence <- gsub("\\.\\.\\.$", "", sequence)  # Remove trailing dots
        extracted <- c(extracted, sequence)
      }
    }
  }

  # Write the extracted sequences to the output file
  writeLines(extracted, output_file)

  cat("Sequences extracted and written to:", output_file, "\n")
}

# Function to format blastclust output
clust2Table <- function(clust, blast) {
  clust_out <- read_tsv(file = clust, col_names = F)
  blast_out <- read_tsv(file = blast, col_names = T)
  ## Count the number of accession numbers in a cluster
  # Counting number of spaces between acc. no. +1
  clust_out$NumAccs <- map(.x = clust_out$X1, function(x) {
    (str_count(string = x, pattern = " ") + 1)
  })
  ## Create empty vectors to store information
  empty_vec <- c("ClusterID")
  empty_vec2 <- c("RowNum")
  ## Adding empty vectors to dataframe
  clust_out[, empty_vec] <- NA
  clust_out[, empty_vec2] <- NA
  ## Counting number of rows to add to the RowNum column -- used for creating cluster name
  rows <- as.numeric(rownames(clust_out)) %>%
    str_pad(width = 4, pad = 0)
  ## Add row number info to dataframe
  clust_out[, "RowNum"] <- rows
  # Name columns
  colnam <- list("AccNum", "NoAccs", "ClusterID", "NumRows")
  ## Create cluster name from row number (num of clusters) and num of accessions
  clust_out <- clust_out %>%
    `colnames<-`(colnam) %>%
    mutate(ClusterID = paste0(NumRows, ".", NoAccs))

  # Store data frame column as vector
  myvar <- c("AccNum", "ClusterID")
  ## Add only the 2 colummns wanted to a new varible
  clusters <- clust_out[myvar]
  # Initialize empty data frame
  new_clust <- data.frame(ClusterID = character(0),
                          AccNum = character(0),
                          stringsAsFactors = F)

  ## Assigning each sseqid to a ClusterID
  # Iterate over dataframe
  for (i in 1:nrow(clusters)) {
    # Extract the cluster name for this row
    cname <- clusters$ClusterID[i]
    # Split accession numbers by a space
    # vals is a vecto of all accession numbers in a row
    vals <- clusters$AccNum[i] %>%
      strsplit(split = " ") %>%
      unlist()

    # Iterate over each element in vals
    for (v in vals) {
      # add each accession num w/ corresponding cluster ID to new df
      new_clust[nrow(new_clust) + 1, ] <- c(cname, v)
    }
  }
  # blast_out$AccNum <- gsub("^>", "", blast_out$AccNum)
  # blast_out$AccNum <- paste0(">", blast_out$AccNum)
  blast_clustnames <- merge(blast_out, new_clust, by = "AccNum")

  for (i in 1:nrow(blast_clustnames)) {
    blast_clustnames$ClusterID[i] <- str_c(blast_clustnames$ClusterID[i])
  }

  first_prot <- as.data.frame(word(clust_out$AccNum), word(clust_out$ClusterID))

  ## write the new file as a TSV file
  newarg <- gsub(".bclust.L[0-9][0-9]S[0-9][0-9].tsv", "", clust)
  # accnum + clusterID
  write_tsv(new_clust, file = paste0(newarg, ".clustIDs"), append = F)
  # first protein from every cluster
  write_tsv(first_prot, file = paste0(newarg, ".clust_reps"),
            col_names = F, append = F)
  # cleaned up blast file + clusterIDs
  write_tsv(blast_clustnames, file = paste0(newarg, ".cln.clust.tsv"),
            col_names = T, append = F)
}

# Function to run InterProScan
runIPRScan2 <- function(query_file, prefix, outdir) {

  # Start InterProScan run
  cat("\n######################\n")
  cat("BEGIN INTERPROSCAN RUN\n")
  cat("Input file:", query_file, "\n")
  cat("######################\n")

  # Output file path
  outfile <- file.path(outdir, paste0(prefix, ".iprscan"))

  # Process the input query file with InterProScan
  cat("Now processing", query_file, "\n")

  # Run InterProScan command
  # Construct the command

  # get the path to the interproscan.sh script from the environment
  # variable INTERPROSCAN_CMD, or assume it's on the path if unspecified
  iprscan_cmd <- Sys.getenv("INTERPROSCAN_CMD", unset="interproscan.sh")

  command <- paste(
    iprscan_cmd, "-i",
    shQuote(query_file),
    "-b", shQuote(outfile),
    "-f TSV --cpu", Sys.getenv("INTERPROSCAN_CPUS", "4"),
    "--appl Pfam,MobiDBlite,Phobius,Coils,SignalP_GRAM_POSITIVE,",
    "SignalP_GRAM_NEGATIVE,Hamap,Gene3D,SignalP_EUK"
  )

  # Run the command
  system(command)

  cat("##################\n")
  cat("END OF IPRSCAN RUN\n")
  cat("##################\n")
}

ipr2Linear <- function(ipr, acc2info, prefix) {
  # read in iprscan results
  # duplicate rows in iprscan file
  ipr_in <- read_tsv(ipr, col_names = TRUE) %>%
    mutate(DB.ID = gsub("G3DSA:", "", DB.ID))

  acc2info_out <- fread(input = acc2info, sep = "\t", header = T, fill = T) %>%
    mutate(FullAccNum = gsub("\\|", "", FullAccNum)) %>%
    mutate(FullAccNum = gsub(".*[a-z]", "", FullAccNum))

  # merge ipr file with acc2info file
  ipr_in <- ipr_in %>%
    # remove version number and any other suffices
    mutate(AccNum.noV = gsub("\\.[0-9].*", "", AccNum))

  ipr_tax <- left_join(ipr_in, acc2info_out, by = "AccNum")

  # read in lineage map
  lineage_map <- fread(
    system.file("common_data", "lineage_lookup.txt", package = "MolEvolvR", mustWork = TRUE),
    header = TRUE, fill = TRUE)

  # merge ipr+info w/ lineage
  # both tables have a species column, but only
  # the lineage_map (y) species column is kept
  ipr_tax <- ipr_tax %>%
    mutate(TaxID = as.numeric(TaxID))

  ipr_lin <- left_join(ipr_tax, lineage_map, by = "TaxID") |>
    mutate(Species = Species.y) %>%
    select(-Species.x, -Species.y)

  # add lookup table to iprscan file
  lookup_tbl <- fread(
    system.file("common_data", "cln_lookup_tbl.tsv", package = "MolEvolvR", mustWork = TRUE),
                      sep = "\t", header = TRUE, fill = TRUE) %>%
    distinct()
  if ("AccNum.x" %in% names(ipr_lin)) {
    ipr_lin <- ipr_lin %>%
      # deselect the AccNum.y from the lineage table (y) and set
      # the AccNum.x (x) from the ipr/acc2info tables to simply 'AccNum'
      mutate(AccNum = AccNum.x) %>%
      select(-AccNum.x, -AccNum.y)
  }
  # run add_name f(x) on ipr+lineage dataframe
  ipr_lin <- ipr_lin %>%
    addName() %>%
    mutate(Name = gsub("^_", "", Name))

  # add domarch info to iprscan + lineage df, only keep what's in x
  ipr_cln <- left_join(ipr_lin, lookup_tbl, by = "DB.ID")

  # populate empty description/short name columns
  for (i in 1:nrow(ipr_cln)) {
    if ((is.na(ipr_cln$ShortName[i]) || ipr_cln$ShortName[i] == "") &&
        (is.na(ipr_cln$SignDesc[i]) || ipr_cln$SignDesc[i] == "-")) {

      ipr_cln$SignDesc[i] <- ipr_cln$IPRDesc[i]
      if (length(ipr_cln$LookupTblDesc[i]) != 0) {
        ipr_cln$ShortName[i] <- ipr_cln$LookupTblDesc[i]
      }
    }
    if (is.na(ipr_cln$ShortName[i]) || ipr_cln$ShortName[i] == "") {
      ipr_cln$ShortName[i] <- ipr_cln$SignDesc[i]
    }
    if (is.na(ipr_cln$SignDesc[i]) || ipr_cln$SignDesc[i] == "-") {
      ipr_cln$SignDesc[i] <- ipr_cln$ShortName[i]
    }
  }
  # rename unclear/duplicated columns
  names(ipr_cln)[names(ipr_cln) == "Description.x"] <- "ProteinName"
  names(ipr_cln)[names(ipr_cln) == "Description.y"] <- "LookupTblDesc"
  # deselect the AccNum.noV from the lineage table (y) and set
  # the AccNum.noV.x (x) from the ipr/acc2info tables to simply 'AccNum.noV'
  ipr_cln <- ipr_cln |> dplyr::select(-AccNum.noV.y) |>
    dplyr::mutate(AccNum.noV = AccNum.noV.x) |>
    dplyr::select(-AccNum.noV.x)
  # create label column to use in ipr2viz
  ipr_cln <- ipr_cln %>%
    mutate(Label = strtrim(ShortName, 30)) %>%
    mutate(Label = gsub(", .*", "", Label)) %>%
    mutate(Label = gsub("C-terminal region of a signal",
                        "C-term signal peptide", Label)) %>%
    mutate(Label = gsub("N-terminal region of a signal",
                        "N-term signal peptide", Label)) %>%
    mutate(Label = gsub("Twin arginine translocation \\(T",
                        "Tat signal profile", Label)) %>%
    mutate(Label = gsub("GLUTATHIONE HYDROLASE PROENZYM",
                        "GLUTATHIONE HYDROLASE PROENZYME", Label)) %>%
    mutate(Label = gsub("N-terminal nucleophile aminohy",
                        "Ntn hydrolases", Label)) %>%
    mutate(Label = gsub("Region of a membrane-bound pro",
                        "cytoplasmic reg of mem-bound prot",
                        Label)) %>%
    mutate(ShortName = gsub("Region of a membrane-bound protein predicted to be
    outside the membrane, in the cytoplasm.",
                            "cytoplasmic reg of mem-bound prot", ShortName))
    outdir <- dirname(ipr)
  # write results to file
  write_tsv(ipr_cln, file.path(paste0(outdir, "/",paste0(prefix,".iprscan_cln.tsv"))))
}

ipr2DomArch <- function(infile_ipr, prefix,
                   analysis = c(
                     "Pfam", "SMART", "Phobius",
                     "Gene3D", "TMHMM", "SignalP_GRAM_POSITIVE",
                     "SUPERFAMILY", "MobiDBLite", "TIGRFAM", "PANTHER", "Coils"
                   )) {
  # read in cleaned up iprscan results
  ipr_in <- read_tsv(infile_ipr, col_names = T, col_types = ipr_cln_cols)

  # split dataframe into unique proteins
  x <- split(x = ipr_in, f = ipr_in$AccNum)

  # plan(strategy = "multicore", .skip = T)

  # within each data.table
  domarch <- map(x, function(y) {
    # domarch <- future_map(x, function(y) {
    acc_row <- data.frame(AccNum = y$AccNum[1], stringsAsFactors = F)
    DAs <- data.frame(matrix(nrow = 1, ncol = length(analysis)))
    DA <- y %>%
      group_by(Analysis) %>%
      arrange(StartLoc)
    i <- 1
    for (a in analysis) {
      a_da <- DA %>% filter(Analysis == a)
      if (a == "SignalP_EUK" || a == "SignalP_GRAM_NEGATIVE" ||
          a == "SignalP_GRAM_POSITIVE") {
        var_shortname <- "DB.ID"
      } else {
        var_shortname <- "ShortName"
      }
      var_shortname_sym <- sym(var_shortname)
      a_da <- a_da %>%
        ungroup() %>%
        select({{ var_shortname_sym }}) %>%
        filter(!is.na({{ var_shortname_sym }})) %>%
        filter(!is.null({{ var_shortname_sym }})) %>%
        pull(var_shortname) %>%
        paste(collapse = "+")
      DAs[1, i] <- a_da
      i <- (i + 1)
    }

    colnames(DAs) <- paste("DomArch", analysis, sep = ".")
    return(cbind(acc_row, DAs))
  })

  # select relevant rows from ipr input to add to domarch
  ipr_select <- ipr_in %>%
    select(Name, AccNum, Species, TaxID, Lineage, Lineage_long_na,
           Lineage_long, Lineage_med, Lineage_short, ProteinName,
           SourceDB, Completeness, AccNum.noV) %>%
    distinct()

  # combine domarchs to one data frame, merge w/ acc2info
  domarch2 <- do.call(rbind.data.frame, domarch)

  domarch_lins <- domarch2 %>%
    merge(ipr_select, by = "AccNum", all.x = T)

  # save domarch_lins file
  write_tsv(domarch_lins,
            file = paste0(prefix, "/", prefix, ".ipr_domarch.tsv"),
            append = F, na = "NA"
  )

  # return domarch2 dataframe to append to blast results if given
  return(domarch2)
}

## function to add results from ipr2DomArch to blast results
appendIPR <- function(ipr_da, blast, prefix) {
  # ! an 'AccNum' or 'AccNum.noV' column is required in blast table for joining !
  blast_out <- read_tsv(blast, col_names = T)
  if ("AccNum.noV" %in% colnames(blast_out)) {
    ipr_da <- read_tsv(paste0(prefix, "/", prefix, ".ipr_domarch.tsv"), col_names = T)
    blast_ipr <- merge(blast_out, ipr_da, by = "AccNum.noV", all.x = T)
  } else {
    blast_ipr <- merge(blast_out, ipr_da, by = "AccNum", all.x = T)
  }

  write_tsv(blast_ipr, file = paste0(prefix, "/", prefix, ".full_analysis.tsv"), na = "NA")
}

# Web BLAST output
web_blast_colnames <- c("Query", "AccNum",
                        "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                        "QStart", "QEnd", "SStart", "SEnd",
                        "EValue", "BitScore", "PcPosOrig",
                        "QSFrames") # specific to "blastx"


# BLAST Command line
cl_blast_colnames <- c("Query", "SAccNum", "AccNum",
                       "SAllSeqID", "STitle", "Species", "TaxID",
                       "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                       "QStart", "QEnd", "QLength",
                       "SStart", "SEnd", "SLength",
                       "EValue", "BitScore", "PcPosOrig",
                       "PcPositive", "ClusterID") # post-cleanup

# IPRSCAN (web+command-line)
ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                  "Status", "RunDate", "IPRAcc", "IPRDesc")

# RPSBLAST
rps_colnames <- c("AccNum", "DBID", "DBSeqID",
                  "PcIdentity", "PcPosOrig", # Ppos missing
                  "AlnLength", "Mismatch",
                  # Q here is Subject; S here is the matching domain; rename!
                  "SStart", "SEnd", "DStart", "DEnd",
                  "EValue", "BitScore", "TaxID") # TaxID missing (NA); remove?

# IPG
ipg_colnames <- c("IPG.ID", "Source", "NucAccNum",
                  "NucStart", "NucStop", "Strand",
                  "AccNum", "ProtDesc",
                  "Species", "SppStrain", "AssemblyID")
# Final ColNames
combo_colnames <- c("Query", "AccNum", "Species", "TaxID", "Lineage",
                    "PcPositive", "ClusterID",
                    # "Leaf", # MISSING (useful for all dataviz)
                    # "AssemblyID", "GeneName", "ProtDesc", # MISSING NOW!?!
                    "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
                    "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP")


## ############ ##
## COLUMN names ##
## ############ ##
## BLAST
############
## Web-BLAST
############
## Downloaded as HIT-TABLE csv
# BLASTP and related protein BLASTs
web_blastp_hit_colnames <- c(
  "Query", "AccNum",
  "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
  "QStart", "QEnd", "SStart", "SEnd",
  "EValue", "BitScore", "PcPosOrig"
)
# BLASTX
web_blastx_colnames <- c(
  "Query", "AccNum",
  "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
  "QStart", "QEnd", "SStart", "SEnd",
  "EValue", "BitScore", "PcPosOrig",
  "QSFrames"
) # specific to "blastx"

## Downloaded as Descriptions csv
# BLASTP and related protein BLASTs
web_blastp_desc_colnames <- c(
  "Description", "Species", "CommonName", "TaxID",
  "BitScore", "TotalScore",
  "PcQCover", "EValue", "PcIdentity",
  "SLen", "AccNum"
)

#####################
## Command line BLAST
#####################

# pre-cleanup
cl_blast_colnames <- c(
  "Query", "SAccNum", "AccNum",
  "SAllSeqID", "STitle", "Species", "TaxID",
  "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
  "QStart", "QEnd", "QLength",
  "SStart", "SEnd", "SLength",
  "EValue", "BitScore", "PcPosOrig"
)

# post-cleanup
cl_blast_postcln_cols <- c(
  "Query", "AccNum",
  "STitle", "Species", "TaxID", "Lineage", "Lineage_long",
  "Lineage_long_na", "Lineage_med", "Lineage_short",
  "PcPositive", "PcIdentity", "AlnLength",
  "SAccNum", "SAllSeqID",
  "Mismatch", "GapOpen",
  "QStart", "QEnd", "QLength",
  "SStart", "SEnd", "SLength",
  "EValue", "BitScore", "PcPosOrig", "QueryName"
)

##########
## IPRSCAN
##########
# ipr_colnames_orig <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
#                        "DB_ID", "SignDesc", "StartLoc", "StopLoc", "Score",
#                        "Status", "RunDate", "IPRAcc", "IPRDesc")

ipr_colnames <- c(
  "AccNum", "SeqMD5Digest", "SLength", "Analysis",
  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
  "Status", "RunDate", "IPRAcc", "IPRDesc"
)

# post cleanup
##########################
## NEED TO BE REORDERED ##
##########################
ipr_cln_colnames <- c(
  "DB.ID", "TaxID", "AccNum.noV", "AccNum",
  "SeqMD5Digest", "SLength", "Analysis", "SignDesc",
  "StartLoc", "StopLoc", "Score", "Status", "RunDate",
  "IPRAcc", "IPRDesc", "FullAccNum", "ProteinName",
  "Length", "SourceDB", "Completeness", "Lineage",
  "Species", "Name", "ShortName", "LookupTblDesc",
  "ID", "Label"
)

###########
## RPSBLAST
###########
rps_colnames <- c(
  "AccNum", "DB.ID", "DBSeqID",
  "PcIdentity.Dom", "PcPosOrig.Dom", # "PcPos.Dom", # Ppos missing
  "AlnLength", "Mismatch",
  "SStart", "SEnd", "DStart", "DEnd",
  "EValue", "BitScore", "TaxID"
) # TaxID missing (NA); remove?

#######################
## IPG and Lineage maps
#######################
ipg_colnames <- c(
  "IPG.ID", "Source", "NucAccNum",
  "NucStart", "NucStop", "Strand",
  "AccNum", "Description",
  "Species", "Spp.Strain", "AssemblyID"
)

##################
## Assembly files
## Genbank, Refseq
##################
assembly_colnames <- c(
  "AssemblyID",
  "bioproject", "biosample", "wgs_master", # not used
  "RefseqCategory", "TaxID", "Spp.TaxID",
  "Species", "Spp.Strain",
  "isolate", "version_status", # not used
  "assembly_level", "release_type", # not used
  "GenomeStatus",
  "seq_rel_date", "asm_name", "submitter", # not used
  "AssemblyID.GBRS",
  "paired_asm_comp", "ftp_path", # not used
  "excluded_from_refseq", "relation_to_type_material"
) # not used
assembly_sub_colnames <- c(
  "TaxID", "Spp.TaxID", "Species", "Spp.Strain",
  "RefseqCategory", "GenomeStatus",
  "AssemblyID", "AssemblyID.GBRS"
)

#################
## Lookup tables
## in common_data
#################
lineage_lookup_colnames <- c("TaxID", "Species", "Lineage_long",
                             "Lineage_long_na", "Lineage_med",
                             "Lineage_short", "Lineage")
domarch_lookup_colnames <- c("DB.ID", "ShortName", "Description", "ID")
# !! SC and LS will fix other piecemeal files based on these

######################
## FINAL UPLOADED DATA
######################
## Combined data frame that is loaded on to the webapp
combo_colnames <- c(
  "Query", "UID", "AccNum", "Species", "TaxID", "Lineage",
  "PcPositive", "ClusterID", "QueryName",
  # "AssemblyID", "GeneName", "Description", # MISSING NOW!?!
  "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
  "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP",
  "DomArch.SMART", "DomArch.TIGR"
)


################
## read tsv colnames
################
lookup_table_cols <- cols(
  DB.ID = col_character(),
  ShortName = col_character(),
  Description = col_character(),
  ID = col_character()
)

iprscan_cols <- cols(
  .default = col_character(),
  TaxID = col_double(),
  SLength = col_double(),
  SignDesc = col_character(),
  StartLoc = col_double(),
  StopLoc = col_double(),
  Score = col_double(),
  Status = col_logical(),
  IPRAcc = col_character(),
  IPRDesc = col_character(),
  Length = col_double(),
  ShortName = col_character(),
  LookupTblDesc = col_character(),
  ID = col_character(),
  Label = col_character()
)

ipr_cln_cols <- cols(
  .default = col_character(),
  TaxID = col_double(),
  SLength = col_double(),
  StartLoc = col_double(),
  StopLoc = col_double(),
  Score = col_double(),
  Status = col_logical(),
  IPRAcc = col_logical(),
  IPRDesc = col_logical(),
  Length = col_double(),
  ID = col_logical()
)

lineage_map_cols <- c(
  "double",
  "character",
  "character", "character", "character", "character", "character"
)
