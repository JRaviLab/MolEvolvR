## Pre-requisites to generate MSA and Phylogenetic Tree
## Includes the following functions:
## convertAlignment2FA, convert2TitleCase, addLeaves2Alignment
## generateAllAlignments2FA
## convertAlignment2TSV??, convertAccNum2FA??
## Created from addLeaves2Alignment.R, convertAlignment2FA.R, all_aln2fa.R
## Modified: Dec 24, 2019 | Jan 2021
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

api_key <- Sys.getenv("ENTREZ_API_KEY", unset = "YOUR_KEY_HERE")

#################
## Pkgs needed ##
#################
# suppressPackageStartupMessages(library(here))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(rentrez))
# suppressPackageStartupMessages(library(msa))
# suppressPackageStartupMessages(library(furrr))
# suppressPackageStartupMessages(library(future))
# suppressPackageStartupMessages(library(doParallel))
# registerDoParallel(cores = detectCores() - 1)
# # library(seqRFLP)
# conflicted::conflict_prefer("filter", "dplyr")
# conflicted::conflict_prefer("strsplit", "base")

##############################
## Pre-requisite functions ##
##############################

## Function to convert to 'Title Case'
#' Changing case to 'Title Case'
#'
#' @author Andrie, Janani Ravi
#' @description Translate string to Title Case w/ delimitter.
#' @aliases totitle, to_title
#' @usage convert2TitleCase(text, delimitter)
#' @param x Character vector.
#' @param y Delimitter. Default is space (" ").
#' @seealso chartr, toupper, and tolower.
#'
#' @return Character vector with the input strings converted to title case.
#'
#' @export
#'
#' @examples
#' # Convert a single string to title case
#' convert2TitleCase("hello world") # Returns "Hello World"
#'
convert2TitleCase <- function(x, y = " ") {
    s <- strsplit(x, y)[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = y
    )
}

################################
## Function to add leaves to an alignment file
## !! Add DA to leaves?
#' addLeaves2Alignment
#'
#' @author Janani Ravi
#' @keywords alignment, accnum, leaves, lineage, species
#' @description Adding Leaves to an alignment file w/ accessions
#' Genomic Contexts vs Domain Architectures.
#'
#' @param aln_file haracter. Path to file. Input tab-delimited file +
#'  alignment file accnum & alignment.
#'  Default is 'pspa_snf7.aln'
#' @param lin_file Character. Path to file. Protein file with accession +
#' number to lineage mapping.
#' Default is 'pspa.txt'
#' @param reduced Boolean. If TRUE, a reduced data frame will be generated with
#' only one sequence per lineage. Default is FALSE.
#'
#' @importFrom dplyr distinct left_join mutate select
#' @importFrom purrr map
#' @importFrom readr read_file read_tsv
#' @importFrom stringr str_sub
#' @importFrom tidyr replace_na separate
#'
#' @return A data frame containing the enriched alignment data with lineage
#' information.
#'
#' @details The alignment file would need two columns: 1. accession +
#' number and 2. alignment. The protein homolog accession to lineage mapping +
#' file should have
#'
#' @note Please refer to the source code if you have alternate +
#' file formats and/or column names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' addLeaves2Alignment("pspa_snf7.aln", "pspa.txt")
#' }
addLeaves2Alignment <- function(aln_file = "",
    lin_file = "data/rawdata_tsv/all_semiclean.txt", # !! finally change to all_clean.txt!!
    # lin_file="data/rawdata_tsv/PspA.txt",
    reduced = FALSE) {
    ## SAMPLE ARGS
    # aln_file <- "data/rawdata_aln/pspc.gismo.aln"
    # lin_file <- "data/rawdata_tsv/all_semiclean.txt"
    # reduced=F

    # 1. Read aln & lineage master files files w/ read_file/read_tsv
    # 2. paste and collapse files so they can be read w/ tsv
    # 3. If the file has 1 column, separate it
    aln <- read_file(aln_file)
    lin <- read_tsv(lin_file)
    aln <- paste(aln, sep = "\\s+", collapse = "\\t")
    aln <- read_tsv(aln, col_names = F)
    if (length(aln) == 1) {
        colnames(aln) <- "x1"
        aln <- separate(aln,
            col = .data$x1,
            into = c("x1", "x2"),
            sep = "\\s+"
        )
    }

    colnames(aln) <- c("AccNum", "Alignment")

    aln_lin <- left_join(aln, lin, by = "AccNum") %>%
        select(
            .data$AccNum, .data$Alignment,
            .data$Species, .data$Lineage
        )

    # Removes rows with NA
    aln_lin <- aln_lin[complete.cases(aln_lin), ]
    # Removes duplicated rows
    aln_lin <- aln_lin %>% distinct()

    ## !! REVISE or REMOVE??? MSA doesn't work with too many sequences!!
    ## !! FIX ASAP.
    if (reduced) {
        # Removes duplicate lineages
        aln_lin <- aln_lin %>% distinct(.data$Lineage, .keep_all = TRUE)
    }

    temp <- aln_lin %>%
        separate(.data$Lineage,
            into = c("Kingdom", "Phylum"),
            sep = ">", remove = F, ## !! How to deal w/ lineages without a phylum?
            extra = "merge", fill = "right"
        ) %>%
        replace_na(replace = list(Kingdom = "", Phylum = "")) %>%
        separate(.data$Species,
            into = c("Genus", "Spp"),
            sep = " ", remove = F,
            extra = "merge", fill = "left"
        ) %>%
        # 3char from kingdom, 6 char from phylum, 1 char from Genus, 3 char from species
        # kingdomPhylum_GenusSpecies
        mutate(Leaf = paste(
            paste0(
                str_sub(.data$Kingdom,
                    start = 1, end = 1
                ),
                str_sub(.data$Phylum, 1, 6)
            ),
            paste0(
                str_sub(.data$Genus, start = 1, end = 1),
                str_sub(.data$Spp, start = 1, end = 3)
            ),
            # AccNum,
            sep = "_"
        ))
    temp$Leaf <- map(temp$Leaf, convert2TitleCase)
    temp <- temp %>%
        mutate(Leaf_Acc = (paste(.data$Leaf, .data$AccNum, sep = "_")))

    # Combine and run through add leaves
    # 3 columns AccNum Sequence Leaf result
    # Create Leaf_AccNum pasted together
    # 2 columns Leaf_AccNum and Sequence Far left
    leaf_aln <- temp %>%
        select(.data$Leaf_Acc, .data$Alignment)
    return(leaf_aln)
}


#' addName
#'
#' @author Samuel Chen, Janani Ravi
#' @description This function adds a new 'Name' column that is comprised of components from
#' Kingdom, Phylum, Genus, and species, as well as the accession
#'
#' @param data Data to add name column to
#' @param accnum_col Column containing accession numbers
#' @param spec_col Column containing species
#' @param lin_col Column containing lineage
#' @param lin_sep Character separating lineage levels
#' @param out_col Column that contains the new 'Name' derived from Species,
#'  Lineage, and AccNum info
#'
#' @importFrom data.table setnames
#' @importFrom dplyr mutate pull select
#' @importFrom tidyr separate
#' @importFrom rlang sym
#' @importFrom stringi stri_replace_all_regex
#' @importFrom stringr word
#'
#' @return Original data with a 'Name' column
#' @export
#'
#' @examples
#' # Example usage of the addName function
#' data <- data.frame(
#'   AccNum = c("ACC123", "ACC456"),
#'   Species = c("Homo sapiens", "Mus musculus"),
#'   Lineage = c("Eukaryota>Chordata", "Eukaryota>Chordata")
#' )
#' enriched_data <- addName(data)
#' enriched_data
addName <- function(data,
    accnum_col = "AccNum", spec_col = "Species", lin_col = "Lineage",
    lin_sep = ">", out_col = "Name") {
    cols <- c(accnum_col, "Kingdom", "Phylum", "Genus", "Spp")
    split_data <- data %>%
        separate(
            col = lin_col, into = c("Kingdom", "Phylum"),
            sep = lin_sep
        ) %>%
        mutate(
            Kingdom = strtrim(.data$Kingdom, 1),
            Phylum = strtrim(.data$Phylum, 6)
        )
    if (!is.null(spec_col)) {
        split_data <- split_data %>%
            separate(spec_col, into = c("Genus", "Spp"), sep = " ") %>%
            mutate(
                Genus = strtrim(.data$Genus, 1),
                Spp = word(string = .data$Spp, start = 1)
            )
    } else {
        split_data$Genus <- ""
        split_data$Spp <- ""
    }

    split_data <- split_data %>%
        select(all_of(cols))
    split_data[is.na(split_data)] <- ""

    accnum_sym <- sym(accnum_col)

    Leaf <- split_data %>%
        mutate(Leaf = paste0(
          .data$Kingdom, .data$Phylum, "_",
          .data$Genus, .data$Spp, "_",
            {{ accnum_sym }}
        )) %>%
        pull(Leaf) %>%
        stringi::stri_replace_all_regex(pattern = "^_|_$", replacement = "") %>%
        stringi::stri_replace_all_regex(pattern = "_+", replacement = "_")

    # out_col=sym(out_col)
    data <- mutate(data, l = Leaf) %>%
        setnames(old = "l", new = out_col)
    return(data)
}

################################
## Function to convert alignment 'aln' to fasta format for MSA + Tree
#' convertAlignment2FA
#'
#' @author Janani Ravi
#' @keywords alignment, accnum, leaves, lineage, species
#' @description Adding Leaves to an alignment file w/ accessions
#' Genomic Contexts vs Domain Architectures.
#'
#' @param aln_file Character. Path to file. Input tab-delimited file +
#'  alignment file accnum & alignment.
#'  Default is 'pspa_snf7.aln'
#' @param lin_file Character. Path to file. Protein file with accession +
#' number to lineage mapping.
#' Default is 'pspa.txt'
#' @param fa_outpath Character. Path to the written fasta file.
#' Default is 'NULL'
#' @param reduced Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
#' Default is 'FALSE'
#'
#' @importFrom readr write_file
#' @importFrom data.table data.table
#'
#' @details The alignment file would need two columns: 1. accession +
#' number and 2. alignment. The protein homolog accession to lineage mapping +
#' file should have
#' @note Please refer to the source code if you have alternate +
#' file formats and/or column names.
#'
#' @return A character string representing the FASTA formatted sequences.
#' If `fa_outpath` is provided, the FASTA will also be saved to the specified
#' file.
#' @export
#'
#' @examples
#' \dontrun{
#' addLeaves2Alignment("pspa_snf7.aln", "pspa.txt")
#' }
#'
convertAlignment2FA <- function(aln_file = "",
    lin_file = "data/rawdata_tsv/all_semiclean.txt", # !! finally change to all_clean.txt!!
    fa_outpath = "",
    reduced = FALSE) {
    ## SAMPLE ARGS
    # aln_file <- "data/rawdata_aln/pspc.gismo.aln"
    # lin_file <- "data/rawdata_tsv/all_semiclean.txt"
    # reduced=F
    # fa_outpath="data/alns/pspc.fasta"

    ## Add leaves
    aln <- addLeaves2Alignment(
        aln = aln_file,
        lin = lin_file,
        reduced = reduced
    )
    names <- aln$Leaf_Acc
    sequences <- aln$Alignment
    aln <- data.table(names, sequences)

    ## Convert to Fasta
    fasta <- ""
    for (i in 1:length(aln$names)) {
        fasta <- paste0(fasta, ">", aln$names[i], "\n", aln$sequences[i], "\n")
    }

    if (!is.null(fa_outpath)) {
        write_file(fasta, fa_outpath)
    }

    # fasta_file <- dataframe2fas(aln, output_path)
    return(fasta)
}

#' mapAcc2Name
#'
#' @description
#' Default renameFA() replacement function. Maps an accession number to its name
#'
#' @param line The line of a fasta file starting with '>'
#' @param acc2name Data Table containing a column of accession numbers and a
#' name column
#' @param acc_col Name of the column containing Accession numbers
#' @param name_col Name of the column containing the names that the accession
#' numbers
#' are mapped to
#'
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#'
#' @return A character string representing the updated FASTA line, where the
#' accession number is replaced with its corresponding name.
#' @export
#'
#' @examples
#' \dontrun{
#' mapAcc2Name(">P12345 some description", acc2name, "AccNum", "Name")
#' }
mapAcc2Name <- function(line, acc2name, acc_col = "AccNum", name_col = "Name") {
    # change to be the name equivalent to an addNames column
    # Find the first ' '
    end_acc <- str_locate(line, " ")[[1]]

    accnum <- substring(line, 2, end_acc - 1)

    acc_sym <- sym(acc_col)
    name_sym <- sym(name_col)
    name <- as.character((filter(acc2name, {{ acc_sym }} == accnum) %>%
        pull({{ name_sym }}))[1])

    return(paste0(">", name))
}

#' Rename the labels of fasta files
#'
#' @param fa_path Path to fasta file
#' @param outpath Path to write altered fasta file to
#' @param replacement_function Function to apply to lines starting with '>'
#' @param ... Additional arguments to pass to replacement_function
#'
#' @importFrom purrr map
#' @importFrom readr read_lines write_lines
#'
#' @return A character vector of the modified lines in the FASTA file.
#' @export
#'
#' @examples
#' \dontrun{
#' renameFA("path/to/input.fasta",
#' "path/to/output.fasta", mapAcc2Name, acc2name)
#' }
renameFA <- function(fa_path, outpath,
    replacement_function = mapAcc2Name, ...) {
    lines <- read_lines(fa_path)
    res <- map(lines, function(x) {
        if (strtrim(x, 1) == ">") {
            return(replacement_function(line = x, ...))
        } else {
            return(x)
        }
    }) %>% unlist()

    write_lines(res, outpath)

    return(res)
}

################################
## generateAllAlignments2FA
#' generateAllAlignments2FA
#'
#' @description
#' Adding Leaves to an alignment file w/ accessions
#'
#' @keywords alignment, accnum, leaves, lineage, species
#' @description Adding Leaves to all alignment files w/ accessions & DAs?
#'
#' @param aln_path Character. Path to alignment files.
#' Default is 'here("data/rawdata_aln/")'
#' @param fa_outpath Character. Path to file. Master protein file with AccNum &
#' lineages.
#' Default is 'here("data/rawdata_tsv/all_semiclean.txt")'
#' @param lin_file Character. Path to the written fasta file.
#' Default is 'here("data/alns/")'.
#' @param reduced Boolean. If TRUE, the fasta file will contain only one
#' sequence per lineage.
#' Default is 'FALSE'.
#'
#' @importFrom purrr pmap
#' @importFrom stringr str_replace_all
#'
#' @return NULL. The function saves the output FASTA files to the specified
#' directory.
#'
#' @details The alignment files would need two columns separated by spaces:
#' 1. AccNum and 2. alignment. The protein homolog file should have AccNum,
#' Species, Lineages.
#' @note Please refer to the source code if you have alternate + file formats
#' and/or column names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' generateAllAlignments2FA()
#' }
generateAllAlignments2FA <- function(aln_path = here("data/rawdata_aln/"),
    fa_outpath = here("data/alns/"),
    lin_file = here("data/rawdata_tsv/all_semiclean.txt"),
    reduced = F) {
    # library(here)
    # library(tidyverse)
    # aln_path <- here("data/rawdata_aln/")
    # outpath <- here("data/alns/")
    # lin_file <- here("data/rawdata_tsv/all_semiclean.txt")

    aln_filenames <- list.files(path = aln_path, pattern = "*.aln")
    aln_filepaths <- paste0(aln_path, "/", aln_filenames)
    variable <- str_replace_all(basename(aln_filenames),
        pattern = ".aln", replacement = ""
    )

    ## Using purrr::pmap
    aln2fa_args <- list(
        aln_file = aln_filepaths,
        fa_outpath = paste0(fa_outpath, "/", variable, ".fa")
    )
    pmap(
        .l = aln2fa_args, .f = convertAlignment2FA,
        lin_file = lin_file,
        reduced = reduced
    )
}


# accessions <- c("P12345","Q9UHC1","O15530","Q14624","P0DTD1")
# accessions <- rep("ANY95992.1", 201)
#' acc2FA
#'
#' @description
#' converts protein accession numbers to a fasta format. Resulting
#' fasta file is written to the outpath.
#'
#' @author Samuel Chen, Janani Ravi
#' @keywords accnum, fasta
#'
#' @param accessions  Character vector containing protein accession numbers to
#' generate fasta sequences for.
#' Function may not work for vectors of length > 10,000
#' @param outpath [str] Location where fasta file should be written to.
#' @param plan Character string specifying the parallel processing strategy to
#' use with the `future` package. Default is "sequential".
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom future future plan value
#' @importFrom purrr map
#' @importFrom rentrez entrez_fetch
#'
#' @return A logical value indicating whether the retrieval and conversion were
#' successful. Returns `TRUE` if successful and `FALSE` otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' acc2FA(accessions = c("ACU53894.1", "APJ14606.1", "ABK37082.1"),
#' outpath = "my_proteins.fasta")
#' Entrez:accessions <- rep("ANY95992.1", 201) |> acc2FA(outpath = "entrez.fa")
#' EBI:accessions <- c("P12345", "Q9UHC1",
#' "O15530", "Q14624", "P0DTD1") |> acc2FA(outpath = "ebi.fa")
#' }
acc2FA <- function(accessions, outpath, plan = "sequential") {
    # validation
    stopifnot(length(accessions) > 0)

    # vector to partitioned list
    partition <- function(v, groups) {
        # partition data to limit number of accessions per POST
        l <- length(v)

        partitioned <- list()
        for (i in 1:groups)
        {
            partitioned[[i]] <- v[seq.int(i, l, groups)]
        }

        return(partitioned)
    }

    # configure futures library
    plan(strategy = plan, .skip = T)
    # group accession numbers into batches of maximum size 200
    # the request header size allowance is unknown, but after some
    # informal testing (molevol_scripts/tests) 200 accession numbers per POST
    # seems a good compromise between efficiency of data per POST and staying
    # under the maximum POST data limit
    groups <- ifelse(
        length(accessions) / 200 >= 1,
        ceiling(length(accessions) / 200),
        1
    )
    accessions_partitioned <- partition(accessions, groups)

    # open file connection
    sink(outpath)

    # construct futures for each group of accession numbers
    a <- map(1:length(accessions_partitioned), function(x) {
        # limit of 10 POSTs/sec w/ key
        if (x %% 10 == 0) {
            Sys.sleep(1)
        }
        f <- future::future(
            rentrez::entrez_fetch(
                id = accessions_partitioned[[x]],
                db = "protein",
                rettype = "fasta",
                api_key = Sys.getenv("ENTREZ_API_KEY")
            )
        )
    })
    # perform POST&GET ('entrez_fetch') request(s); print output (output is directed to `outpath`)
    for (f in a)
    {
        cat(future::value(f))
    }
    # close file connection
    sink(NULL)
    # validate the result
    result <- tryCatch(
        expr = {
            Biostrings::readAAStringSet(outpath)
            TRUE
        },
        error = function(e) {
            print(e)
            FALSE
        }
    )
    return(result)
}

#' createRepresentativeAccNum
#'
#' @description
#' Function to generate a vector of one Accession number per distinct
#' observation from 'reduced' column
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @param prot_data Data frame containing Accession Numbers
#' @param reduced Column from prot_data from which distinct observations
#' will be generated from.
#' One accession number will be assigned for each of these observations
#' @param accnum_col Column from prot_data that contains Accession Numbers
#'
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#'
#' @return A character vector containing one Accession number per distinct
#' observation from the specified reduced column.
#' @export
#'
#' @examples
#' \dontrun{
#' createRepresentativeAccNum(prot)
#' }
createRepresentativeAccNum <- function(prot_data,
    reduced = "Lineage",
    accnum_col = "AccNum") {
    # Get Unique reduced column and then bind the AccNums back to get one
    # AccNum per reduced column
    reduced_sym <- sym(reduced)
    accnum_sym <- sym(accnum_col)

    distinct_reduced <- prot_data %>%
        pull(reduced) %>%
        unique()

    accessions <- c()

    for (r in distinct_reduced)
    {
        if (is.na(r)) {
            next
        }
        # r <- toString(distinct_reduced[i,reduced])

        AccNum <- prot_data %>%
            filter({{ reduced_sym }} == r) %>%
            pull(accnum_col)

        accessions <- append(accessions, AccNum[1])
    }
    return(accessions)
}

#' alignFasta
#'
#' @description
#' Perform a Multiple Sequence Alignment on a FASTA file.
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @param fasta_file Path to the FASTA file to be aligned
#' @param tool Type of alignment tool to use. One of three options: "Muscle",
#' "ClustalO", or "ClustalW"
#' @param outpath Path to write the resulting alignment to as a FASTA file.
#' If NULL, no file is written
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom msa msaClustalOmega msaMuscle msaClustalW
#'
#' @return aligned fasta sequence as a MsaAAMultipleAlignment object
#' @export
#'
#' @examples
#' \dontrun{
#' aligned_sequences <- alignFasta("my_sequences.fasta",
#' tool = "Muscle", outpath = "aligned_output.fasta")
#' }
alignFasta <- function(fasta_file, tool = "Muscle", outpath = NULL) {
    fasta <- readAAStringSet(fasta_file)

    aligned <- switch(tool,
        "Muscle" = msaMuscle(fasta),
        "ClustalO" = msaClustalOmega(fasta),
        "ClustalW" = msaClustalW(fasta)
    )

    if (typeof(outpath) == "character") {
        writeMSA_AA2FA(aligned, outpath)
    }
    return(aligned)
}

#' writeMSA_AA2FA
#'
#' @description
#' MsaAAMultipleAlignment Objects are generated from calls to msaClustalOmega
#' and msaMuscle from the 'msa' package
#' @author Samuel Chen, Janani Ravi
#'
#' @param alignment MsaAAMultipleAlignment object to be written as a fasta
#' @param outpath Where the resulting FASTA file should be written to
#'
#' @importFrom Biostrings toString unmasked
#' @importFrom readr write_file
#'
#' @return Character string representing the content of the written FASTA file.
#' @export
#'
#' @examples
#' \dontrun{
#' writeMSA_AA2FA("my_sequences.fasta", outpath = "aligned_output.fasta")
#' }
writeMSA_AA2FA <- function(writeMSA_AA2FA, outpath) {
    l <- length(rownames(alignment))
    fasta <- ""
    for (i in 1:l)
    {
        fasta <- paste0(fasta, paste(">", rownames(alignment)[i]), "\n")
        seq <- toString(unmasked(alignment)[[i]])
        fasta <- paste0(fasta, seq, "\n")
    }
    write_file(fasta, outpath)
    return(fasta)
}

#' getAccNumFromFA
#'
#' @param fasta_file Character. The path to the FASTA file from which
#' accession numbers will be extracted.
#'
#' @importFrom stringi stri_extract_all_regex
#'
#' @return A character vector containing the extracted accession numbers.
#' @export
#'
#' @examples
#' \dontrun{
#' getAccNumFromFA("my_sequences.fasta")
#' }
getAccNumFromFA <- function(fasta_file) {
    txt <- read_file(fasta_file)
    accnums <- stringi::stri_extract_all_regex(fasta_file, "(?<=>)[\\w,.]+")[[1]]
    return(accnums)
}


################################
## convertAccNum2FA
#######
## 1 ##
#######
## Accnum2fa
# ref <- c("U15717", "U15718", "U15719", "U15720",
#          "U15721", "U15722", "U15723", "U15724")
# ref_gb <- read.GenBank(ref)
# cbind(attr(ref_gb, "species"), names(ref_gb))
# attr(ref_gb, "description")

#######
## 2 ##
#######
## GenBank 2 fasta file format
# seqinr::gb2fasta(source.file="", destination.file="")

#######
## 3 ##
#######
## Ref: https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
# retrieveseqs <- function(seqnames,acnucdb)
# {
#   myseqs <- list()   # Make a list to store the sequences
#   require("seqinr")  # This function requires the SeqinR R package
#   choosebank(acnucdb)
#   for (i in 1:length(seqnames))
#   {
#     seqname <- seqnames[i]
#     print(paste("Retrieving sequence",seqname,"..."))
#     queryname <- "query2"
#     query <- paste("AC=",seqname,sep="")
#     query(`queryname`,`query`)
#     seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
#     myseqs[[i]] <- seq
#   }
#   closebank()
#   return(myseqs)
# }
# seqnames <- c("Q10572","E3M2K8","Q8WS01","E1FUV2","A8NSK3","Q9VT99")
# seqs <- retrieveseqs(seqnames,"swissprot")

################################
## convertAlignment2TSV
## NEEDS FIXING!
# convertAlignment2TSV <- function(file_path){
#   cfile <- read_delim("data/alignments/pspc.gismo.aln", delim=" ")
#   cfile <- as.data.frame(map(cfile,function(x) gsub("\\s+", "",x)))
#   colnames(cfile) <- c("AccNum", "Alignment")
# }
