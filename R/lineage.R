# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(rentrez))
# suppressPackageStartupMessages(library(future))
# suppressPackageStartupMessages(library(furrr))
# suppressPackageStartupMessages(library(data.table))

#####################################
## Download Assembly Summary Files ##
#####################################
#' Download the combined assembly summaries of genbank and refseq
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @param outpath String of path where the assembly summary file should be
#' written
#' @param keep Character vector containing which columns should be retained and
#' downloaded
#'
#' @importFrom data.table fwrite setnames
#' @importFrom dplyr bind_rows select
#' @importFrom biomartr getKingdomAssemblySummary
#'
#' @return A tab-separated file containing the assembly summary. The function
#' does notreturn any value but writes the output directly to the specified file.
#' @export
#'
#' @examples
#' \dontrun{
#' downloadAssemblySummary(outpath = "assembly_summary.tsv",
#'      keep = c("assembly_accession", "taxid", "organism_name"))
#' }
downloadAssemblySummary <- function(outpath,
    keep = c(
        "assembly_accession", "taxid",
        "species_taxid", "organism_name"
    )) {
    assembly_kingdom_genbank <- getKingdomAssemblySummary("genbank")
    assembly_kingdom_refseq <- getKingdomAssemblySummary("refseq")

    if (keep == "all") {
        assembly_all <- bind_rows(assembly_kingdom_genbank, assembly_kingdom_refseq)
    } else {
        assembly_all <- bind_rows(assembly_kingdom_genbank, assembly_kingdom_refseq) %>%
            select(all_of(keep))
    }

    assembly_all <- assembly_all %>% data.table::setnames(
        old = c(
            "taxid", "refseq_category", "species_taxid", "organism_name",
            "infraspecific_name", "genome_rep"
        ),
        new = c(
            "TaxID", "RefseqCategory", "Parent.TaxID", "Species",
            "Spp.Strain", "GenomeStatus"
        ),
        skip_absent = T
    )

    # dplyr::rename("AssemblyID"="assembly_accession",
    #               "TaxID"="taxid",
    #               "RefseqCategory"="refseq_category",
    #               "Parent.TaxID"="species_taxid",
    #               "Species"="organism_name",
    #               "Spp.Strain"="infraspecific_name",
    #               "GenomeStatus"="genome_rep")

    fwrite(assembly_all, outpath, sep = "\t")
}



###################################
## Map GCA_ID to TaxID & Lineage ##
###################################
#' Function to map GCA_ID to TaxID, and TaxID to Lineage
#'
#' @author Samuel Chen, Janani Ravi
#' @note
#' Currently configured to have at most kingdom and phylum
#'
#'
#' @param prot_data Dataframe containing a column `GCA_ID`
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "downloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "createLineageLookup()" function
#' @param acc_col Character. The name of the column in `prot_data` containing
#' accession numbers. Default is "AccNum".
#'
#' @importFrom dplyr pull
#' @importFrom data.table fread setnames
#'
#' @return A dataframe containing the merged information of GCA_IDs, TaxIDs,
#' and their corresponding lineage up to the phylum level. The dataframe
#' will include information from the input `prot_data` and lineage data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- GCA2Lineage(prot_data = my_prot_data,
#'                        assembly_path = "path/to/assembly_summary.txt",
#'                        lineagelookup_path = "path/to/lineage_lookup.tsv")
#' }
GCA2Lineage <- function(prot_data,
    assembly_path = "/data/research/jravilab/common_data/assembly_summary_genbank.txt",
    lineagelookup_path = "/data/research/jravilab/common_data/lineage_lookup.tsv",
    acc_col = "AccNum") {
    assembly_summary <- fread(assembly_path, sep = "\t")
    assembly_summary <- setnames(assembly_summary, "AssemblyID", "GCA_ID")

    mergedTax <- merge(
        x = prot_data, y = assembly_summary,
        by = "GCA_ID", all.x = T
    )
    accessions <- prot_data %>%
        pull(acc_col) %>%
        unique()
    # Prioritize Complete Genome
    best_rows <- integer(length(accessions))
    for (i in 1:length(accessions))
    {
        # browser()
        acc <- accessions[i]
        acc_inds <- which(mergedTax$Protein == acc)
        if (length(acc_inds) > 1) {
            complete <- acc_inds[which(mergedTax[acc_inds, ]$assembly_level == "Complete Genome")]
            if (length(complete) != 0) {
                best_rows[i] <- complete[1]
            } else {
                best_rows[i] <- acc_inds[1]
            }
        }
    }
    mergedTax <- mergedTax[best_rows, ]

    lineage_map <- fread(lineagelookup_path, sep = "\t")
    lineage_map <- lineage_map[, !"Species"]

    mergedLins <- merge(mergedTax, lineage_map,
        by.x = "TaxID", by.y = "TaxID",
        all.x = T
    )

    return(mergedLins)
}

###################################
## !! @SAM why is this called lins?
###################################
#' addLineage
#'
#' @param df Dataframe containing accession numbers. The dataframe should
#' have a column specified by `acc_col` that contains these accession numbers.
#' @param acc_col Character. The name of the column in `df` containing
#' accession numbers. Default is "AccNum".
#' @param assembly_path String. The path to the assembly summary file generated
#' using the `downloadAssemblySummary()` function.
#' @param lineagelookup_path String. The path to the lineage lookup file (taxid
#' to lineage mapping) generated using the `create_lineage_lookup()` function.
#' @param ipgout_path String. Optional path to save intermediate output files.
#' Default is NULL.
#' @param plan Character. Specifies the execution plan for parallel processing.
#' Default is "multicore".
#'
#' @importFrom dplyr pull
#' @importFrom rlang sym
#'
#' @return A dataframe that combines the original dataframe `df` with lineage
#' information retrieved based on the provided accession numbers.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' enriched_df <- addLineage(df = my_data,
#'                            acc_col = "AccNum",
#'                            assembly_path = "path/to/assembly_summary.txt",
#'                            lineagelookup_path = "path/to/lineage_lookup.tsv")
#' }
addLineage <- function(df, acc_col = "AccNum", assembly_path,
    lineagelookup_path, ipgout_path = NULL, plan = "multicore") {
    acc_sym <- sym(acc_col)
    accessions <- df %>% pull(acc_sym)
    lins <- acc2Lineage(accessions, assembly_path,
        lineagelookup_path, ipgout_path,
        plan = plan
    )

    # Drop a lot of the unimportant columns for now? will make merging much easier
    lins <- lins[, c(
        "Strand", "Start", "Stop", "Nucleotide Accession", "Source",
        "Id", "Strain"
    ) := NULL]
    lins <- unique(lins)

    # dup <- lins %>% group_by(Protein) %>% summarize(count=n()) %>% filter(count > 1) %>%
    #   pull(Protein)

    ## !! @SAM: there is no "Protein" column anymore !!
    merged <- merge(df, lins,
        by.x = acc_sym, by.y = "Protein", all.x = TRUE
    )
    return(merged)
}

#######################################
## Map Protein Accessions to Lineage ##
#######################################
#' acc2Lineage
#'
#' @description
#' Function to map protein accession numbers to lineage
#'
#' @author Samuel Chen, Janani Ravi
#' @description This function combines 'efetchIPG()' and 'IPG2Lineage()' to map a set
#' of protein accessions to their assembly (GCA_ID), tax ID, and lineage.
#'
#' @param accessions Character vector of protein accessions
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "downloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' @param ipgout_path Path to write the results of the efetch run of the accessions
#' on the ipg database. If NULL, the file will not be written. Defaults to NULL
#' @param plan Character. Specifies the execution plan for parallel processing.
#' Default is "multicore".
#'
#' @return A dataframe containing lineage information mapped to the given protein
#' accessions. The dataframe includes relevant columns such as TaxID, GCA_ID,
#' Protein, Protein Name, Species, and Lineage.
#' @export
#'
#' @examples
#' \dontrun{
#' lineage_data <- acc2Lineage(
#'   accessions = c("P12345", "Q67890"),
#'   assembly_path = "path/to/assembly_summary.txt",
#'   lineagelookup_path = "path/to/lineage_lookup.tsv",
#'   ipgout_path = "path/to/output.txt"
#' )
#' }
acc2Lineage <- function(accessions, assembly_path, lineagelookup_path,
    ipgout_path = NULL, plan = "multicore") {
    tmp_ipg <- F

    if (is.null(ipgout_path)) {
        tmp_ipg <- T
        ipgout_path <- tempfile("ipg", fileext = ".txt")
    }
    efetchIPG(accessions, out_path = ipgout_path, plan = plan)

    lins <- IPG2Lineage(accessions, ipgout_path, assembly_path, lineagelookup_path)

    # if(tmp_ipg)
    # {
    #   unlink(tempdir(), recursive=T)
    # }

    # cols <- c("TaxID","GCA_ID", "Protein", "Protein Name", "Species", "Lineage")
    # lins <- unique(lins[,..cols])

    return(lins)
}


#########################################
## Download IPG results for Accessions ##
#########################################
#' efetchIPG
#'
#' @author Samuel Chen, Janani Ravi
#' @description Perform efetch on the ipg database and write the results to out_path
#'
#' @param accessions Character vector containing the accession numbers to query on
#' the ipg database
#' @param out_path Path to write the efetch results to
#' @param plan Character. Specifies the execution plan for parallel processing.
#' Default is "multicore".
#'
#' @importFrom future future plan
#' @importFrom purrr map
#' @importFrom rentrez entrez_fetch
#'
#' @return The function does not return a value but writes the efetch results
#' directly to the specified `out_path`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' efetchIPG(
#'   accessions = c("P12345", "Q67890", "A12345"),
#'   out_path = "path/to/efetch_results.xml"
#' )
#' }
efetchIPG <- function(accessions, out_path, plan = "multicore") {
    if (length(accessions) > 0) {
        partition <- function(v, groups) {
            # Partition data to limit number of queries per second for rentrez fetch:
            # limit of 10/second w/ key
            l <- length(v)

            partitioned <- list()
            for (i in 1:groups)
            {
                partitioned[[i]] <- v[seq.int(i, l, groups)]
            }

            return(partitioned)
        }

        plan(strategy = plan, .skip = T)

        ## ! Note: LS changed it to 600 because she has 5K results and wanted x ≤ 9
        min_groups <- length(accessions) / 600
        groups <- min(max(min_groups, 15), length(accessions))
        partitioned_acc <- partition(accessions, groups)
        sink(out_path)

        a <- map(1:length(partitioned_acc), function(x) {
            # Avoid hitting the rate API limit
            if (plan != "sequential" & x %% 9 == 0) {
                Sys.sleep(1)
            }
            f <- future({
                entrez_fetch(
                    id = partitioned_acc[[x]],
                    db = "ipg",
                    rettype = "xml", # parsed=T,
                    api_key = "YOUR_KEY_HERE"
                )
            })
        })

        for (f in a)
        {
            cat(value(f))
        }
        sink(NULL)
    }
}

#########################################
## Maps IPG results to TaxID + Lineage ##
#########################################
#' IPG2Lineage
#'
#' @author Samuel Chen, Janani Ravi
#' @description Takes the resulting file of an efetch run on the ipg database and
#' append lineage, and taxid columns
#'
#' @param accessions Character vector of protein accessions
#' @param ipg_file Path to the file containing results of an efetch run on the
#' ipg database. The protein accession in 'accessions' should be contained in this
#' file
#' @param refseq_assembly_path String. Path to the RefSeq assembly summary file.
#' @param genbank_assembly_path String. Path to the GenBank assembly summary file.
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "createLineageLookup()" function
#'
#' @importFrom data.table fread setnames
#'
#' @return A data table containing protein accessions along with their
#' corresponding TaxIDs and lineage information.
#' @export
#'
#' @examples
#' \dontrun{
#' lins <- IPG2Lineage(
#'   accessions = c("P12345", "Q67890"),
#'   ipg_file = "path/to/ipg_results.txt",
#'   refseq_assembly_path = "path/to/refseq_assembly_summary.txt",
#'   genbank_assembly_path = "path/to/genbank_assembly_summary.txt",
#'   lineagelookup_path = "path/to/lineage_lookup.tsv"
#' )
#' }
IPG2Lineage <- function(accessions, ipg_file,
    refseq_assembly_path, genbank_assembly_path,
    lineagelookup_path) {
    ipg_dt <- fread(ipg_file, sep = "\t", fill = T)

    accessions <- unique(accessions)
    ipg_dt <- ipg_dt[Protein %in% accessions]

    ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")

    # Call GCA2Lins with different assembly_paths depending on refseq or not
    # Select for Refseq rows over other DB rows
    refseq_rows <- integer(length(accessions))
    genbank_rows <- integer(length(accessions))
    for (i in 1:length(accessions))
    {
        # browser()
        acc <- accessions[i]
        acc_inds <- which(mergedTax$Protein == acc)
        if (length(acc_inds) != 0) {
            # refseq inds take precedence
            refseq_inds <- acc_inds[which(mergedTax[acc_inds, ]$Source == "RefSeq")]
            if (length(refseq_inds) != 0) {
                # Take the first first row of the refseq (smallest index)
                refseq_rows[i] <- refseq_inds[1]
            } else {
                # take the first row of whatever is left?
                genbank_rows[i] <- acc_inds[1]
            }
        }
    }

    # Empty values be gone
    refseq_rows <- refseq_rows[which(refseq_rows != 0)]
    genbank_rows <- genbank_rows[which(genbank_rows != 0)]

    # Call GCA2Lineages using refseq
    ### Possible to run these in parallel if it takes a while
    if (length(refseq_rows) != 0) {
        refseq_ipg_dt <- ipg_dt[refseq_rows, ]
        refseq_lins <- GCA2Lineage(refseq_ipg_dt,
            assembly_path = refseq_assembly_path,
            lineagelookup_path
        )
    }
    if (length(genbank_rows) != 0) {
        genbank_ipg_dt <- ipg_dt[genbank_rows, ]
        genbank_lins <- GCA2Lineage(gca_ipg_dt,
            assembly_path = genbank_assembly_path,
            lineagelookup_path
        )
    }


    lins <- GCA2Lineage(prot_data = ipg_dt, assembly_path, lineagelookup_path)
    lins <- lins[!is.na(Lineage)] %>% unique()

    return(lins)
}


#########################################
## !! @SAM: Add TaxID based on AccNum? ##
#########################################
#' addTaxID
#'
#' @param data A data frame or data table containing protein accession numbers.
#' @param acc_col A string specifying the column name in `data` that contains
#' the accession numbers. Defaults to "AccNum".
#' @param version A logical indicating whether to remove the last two characters
#' from the accession numbers for TaxID retrieval. Defaults to TRUE.
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table merge.data.table
#'
#' @return A data table that includes the original data along with a new column
#' containing the corresponding TaxIDs.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a sample data table with accession numbers
#' sample_data <- data.table(AccNum = c("ABC123.1", "XYZ456.1", "LMN789.2"))
#' enriched_data <- addTaxID(sample_data, acc_col = "AccNum", version = TRUE)
#' enriched_data
#' }
addTaxID <- function(data, acc_col = "AccNum", version = T) {
    if (!is.data.table(data)) {
        data <- as.data.table(data)
    }

    accessions <- data[[acc_col]]

    if (version) {
        data <- data[, AccNum.noV := substr(data[[acc_col]],
            start = 0, stop = nchar(data[[acc_col]]) - 2
        )]
        acc_col <- "AccNum.noV"
    }

    out_path <- tempdir()
    tax <- proteinAcc2TaxID(accessions, "TEMPTAX", out_path, return_dt = TRUE)

    data <- merge.data.table(data, tax,
        by.x = acc_col, by.y = "AccNum.noV", all.x = T
    )
    return(data)
}

##################################
## Maps Protein AccNum to TaxID ##
##################################
#' proteinAcc2TaxID
#'
#' @param accnums A character vector of protein accession numbers to be mapped
#' to TaxIDs.
#' @param suffix A string suffix used to name the output file generated by the
#' script.
#' @param out_path A string specifying the directory where the output file will
#' be saved.
#' @param return_dt A logical indicating whether to return the result as a data
#' table. Defaults to FALSE. If TRUE, the output file is read into a data table
#' and returned.
#'
#' @importFrom data.table fread
#'
#' @return If `return_dt` is TRUE, a data table containing the mapping of protein
#' accession numbers to TaxIDs. If FALSE, the function returns NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example accession numbers
#' accessions <- c("ABC123", "XYZ456", "LMN789")
#' tax_data <- proteinAcc2TaxID(accessions, suffix = "example",
#' out_path = "/path/to/output", return_dt = TRUE)
#' tax_data
#' }
proteinAcc2TaxID <- function(accnums, suffix, out_path, return_dt = FALSE) {
    # Write accnums to a file
    acc_file <- tempfile()
    write(paste(accnums, collapse = "\n"), acc_file)
    script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh"
    call <- paste(script, acc_file, suffix, out_path)
    system(call, wait = TRUE)
    if (return_dt) {
        out_file <- paste0(out_path, "/", suffix, ".acc2info.tsv")
        dt <- fread(out_file, sep = "\t", fill = T)
        return(dt)
    }
}


#######################################
## OLD: Maps Protein AccNum to TaxID ##
#######################################
#' proteinAcc2TaxID_old
#'
#' @author Samuel Chen, Janani Ravi
#' @description Perform elink to go from protein database to taxonomy database
#' and write the resulting file of taxid and lineage to out_path
#'
#' @param accessions A character vector containing the accession numbers to query
#' in the protein database.
#' @param out_path A string specifying the path where the results of the query
#' will be written. If set to NULL, a temporary directory will be used.
#' @param plan A character string that specifies the execution plan for parallel
#' processing. The default is "multicore".
#'
#' @importFrom future plan
#' @importFrom purrr map
#'
#' @return This function does not return a value. It writes the results to the
#'         specified output path.
#' @export
#'
#' @examples
#' \dontrun{
#' accessions <- c("ABC123", "XYZ456", "LMN789")
#' proteinAcc2TaxID_old(accessions, out_path = "/path/to/output")
#' }
proteinAcc2TaxID_old <- function(accessions, out_path, plan = "multicore") {
    if (length(accessions) > 0) {
        partition <- function(v, groups) {
            # Partition data to limit number of queries per second for rentrez fetch:
            # limit of 10/second w/ key
            l <- length(v)

            partitioned <- list()
            for (i in 1:groups)
            {
                partitioned[[i]] <- v[seq.int(i, l, groups)]
            }

            return(partitioned)
        }

        plan(strategy = plan, .skip = T)

        ## ! Note: LS changed it to 600 because she has 5K results and wanted x to be ≤ 9
        min_groups <- length(accessions) / 600
        groups <- min(max(min_groups, 15), length(accessions))
        partitioned_acc <- partition(accessions, groups)

        out_path <- tempdir()

        a <- map(1:length(partitioned_acc), function(x) {
            # Avoid hitting the rate API limit
            if (plan != "sequential" & x %% 9 == 0) {
                Sys.sleep(1)
            }
            print(x)
            script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh"
            # script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/proteinAcc2TaxID.sh"

            # accnum_in <- paste(partitioned_acc[[x]], collapse=",")
            accnum_in <- tempfile()
            write(paste(partitioned_acc[[x]], collapse = ","), accnum_in)

            system(call, wait = F)
            # system(paste(script, accnum_in), wait=TRUE)


            # f <- future({
            #   el=entrez_link(dbfrom="protein", id=partitioned_acc[[x]],
            #                    db="taxonomy",
            #                    by_id=FALSE,
            #                    api_key="YOUR_KEY_HERE")
            #   entrez_fetch(db="taxonomy",
            #                id=el, rettype="taxid",
            #                api_key="YOUR_KEY_HERE")
            #   # Calling Janani's shell script would be easier
            # })
        })

        # for( f in a)
        # {
        # cat(value(f))
        # }
        # sink(NULL)
    }
}
