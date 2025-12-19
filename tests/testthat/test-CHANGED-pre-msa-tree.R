context("CHANGED-test-pre-msa-tree")
test_that("CHANGED-test-pre-msa-tree", {
    # Test with standard input
    expect_equal(convert2TitleCase("hello world"), "Hello World")
    expect_equal(convert2TitleCase("this is a test"), "This Is A Test")
    
    # Test with single word
    expect_equal(convert2TitleCase("R"), "R")
    expect_equal(convert2TitleCase("hello"), "Hello")
    
    # Test with punctuation
    expect_equal(convert2TitleCase("capitalize every word."), 
                 "Capitalize Every Word.")
    expect_equal(convert2TitleCase("test, case, example"), 
                 "Test, Case, Example")
    
    # Test with empty string
    expect_equal(convert2TitleCase(""), "")
    
    # Test with multiple spaces
    expect_equal(convert2TitleCase("  multiple   spaces  "), 
                 "  Multiple   Spaces ")
    
    # Create mock alignment data
    mock_alignment_data <- c(
        "ACU53894.1\tMKTIIALSYI",
        "APJ14606.1\tGSLKVSISY",
        "ABK37082.1\tMKTLVAILA"
    )
    
    # Write the mock alignment data to a .aln file
    aln_file <- file.path(tempdir(), "mock_alignment.aln")
    writeLines(mock_alignment_data, aln_file)
    
    # Create mock lineage data
    mock_lineage_data <- data.frame(
        AccNum = c("ACU53894.1", "APJ14606.1", "ABK37082.1"),
        Species = c("Species1", "Species2", "Species3"),
        Lineage = c("Kingdom>Phylum1", "Kingdom>Phylum2", "Kingdom>Phylum3"),
        stringsAsFactors = FALSE
    )
    
    # Write lineage data to a TSV file
    lin_file <- file.path(tempdir(), "mock_lineage.tsv")
    write_tsv(mock_lineage_data, lin_file)
    
    # Now you can test the function
    result <- addLeaves2Alignment(aln_file, lin_file, seq_type = 'AA')
    
    # Check the structure of the result
    expect_is(result, "data.frame")
    
    # Check the columns in the result
    expect_true("Leaf_Acc" %in% colnames(result))
    expect_true("Alignment" %in% colnames(result))
    
    # Check the content of the Leaf_Acc column
    expect_equal(result$Leaf_Acc[1], "KPhylum_NASpe_ACU53894.1") 
    expect_equal(result$Leaf_Acc[2], "KPhylum_NASpe_APJ14606.1")   
    expect_equal(result$Leaf_Acc[3], "KPhylum_NASpe_ABK37082.1") 
    
    # Check alignment sequences correspond to accession numbers
    expect_equal(result$Alignment[1], "MKTIIALSYI")
    expect_equal(result$Alignment[2], "GSLKVSISY")
    expect_equal(result$Alignment[3], "MKTLVAILA")
    
    # Mock data to simulate input
    mock_data <- tibble::tibble(
        AccNum = c("Acc1", "Acc2", "Acc3", "Acc4"),
        Species = c("Homo sapiens", "Escherichia coli", "Mus musculus", NA),
        Lineage = c("Eukaryota>Chordata", "Bacteria>Proteobacteria", 
                    "Eukaryota>Chordata", "Bacteria>Firmicutes")
    )
    result <- addName(mock_data)
    
    # Check that the output contains the new Name column
    expect_true("Name" %in% colnames(result))
    
    # Check if the names are formatted correctly
    expected_names <- c("EChorda_Hsapiens_Acc1", "BProteo_Ecoli_Acc2", 
                        "EChorda_Mmusculus_Acc3", "BFirmic_Acc4")
    expect_equal(result$Name, expected_names)
    
    # Test with NULL species column
    result <- addName(mock_data, spec_col = NULL)
    
    # Check that species is ignored, and names are based on lineage and AccNum
    expected_names <- c("EChorda_Acc1", "BProteo_Acc2",
                        "EChorda_Acc3", "BFirmic_Acc4")
    expect_equal(result$Name, expected_names)
    
    # Create mock data with missing values
    mock_data_with_missing <- tibble::tibble(
        AccNum = c("Acc1", "Acc2", "Acc3"),
        Species = c(NA, "Escherichia coli", NA),
        Lineage = c(NA, "Bacteria>Proteobacteria", "Eukaryota>Chordata")
    )
    
    result <- addName(mock_data_with_missing)
    
    # Check that missing species and lineage data are handled correctly
    expected_names <- c("_Acc1", "BProteo_Ecoli_Acc2", "EChorda_Acc3")
    expect_equal(result$Name, expected_names)
    
    # Test with a different separator (e.g., "/")
    mock_data_custom_sep <- mock_data
    mock_data_custom_sep$Lineage <- gsub(">", "/", mock_data_custom_sep$Lineage)
    
    result <- addName(mock_data_custom_sep, lin_sep = "/")
    
    # Ensure the function correctly parses the lineage with the custom separator
    expected_names <- c("EChorda_Hsapiens_Acc1", 
                        "BProteo_Ecoli_Acc2", 
                        "EChorda_Mmusculus_Acc3", 
                        "BFirmic_Acc4")
    expect_equal(result$Name, expected_names)
    
    # Test with an empty dataframe
    empty_data <- tibble::tibble(AccNum = character(0), 
                                 Species = character(0), 
                                 Lineage = character(0))
    
    result <- addName(empty_data)
    
    # Ensure the result is still a dataframe with the expected structure
    expect_equal(nrow(result), 0)
    expect_true("Name" %in% colnames(result))
    
    # Create mock data that might result in extra underscores
    mock_data_with_gaps <- tibble::tibble(
        AccNum = c("Acc1", "Acc2"),
        Species = c(NA, "Escherichia coli"),
        Lineage = c("Eukaryota", "Bacteria>Proteobacteria")
    )
    
    result <- addName(mock_data_with_gaps)
    
    # Ensure no extra underscores are present in the Name column
    expected_names <- c("E_Acc1", "BProteo_Ecoli_Acc2")
    expect_equal(result$Name, expected_names)
    


    
    # Call the function to convert alignment to FASTA
    result <- convertAlignment2FA(
        aln_file = aln_file,
        lin_file = lin_file,
        fa_outpath = NULL,
        reduced = FALSE
    )
    
    # Expected FASTA format output
    expected_fasta <- ">KPhylum_NASpe_ACU53894.1\nMKTIIALSYI\n>KPhylum_NASpe_APJ14606.1\nGSLKVSISY\n>KPhylum_NASpe_ABK37082.1\nMKTLVAILA\n"
    
    # Check if the result matches the expected FASTA output
    expect_equal(result, expected_fasta)
    
    # Temporary output file path
    fa_outpath <- tempfile(fileext = ".fasta")
    
    # Call the function to convert alignment to FASTA and save it to a file
    convertAlignment2FA(
        aln_file = aln_file,
        lin_file = lin_file,
        fa_outpath = fa_outpath,
        reduced = FALSE
    )
    
    # Read the generated FASTA file
    fasta_content <- read_file(fa_outpath)
    
    # Expected FASTA format output
    expected_fasta <- ">KPhylum_NASpe_ACU53894.1\nMKTIIALSYI\n>KPhylum_NASpe_APJ14606.1\nGSLKVSISY\n>KPhylum_NASpe_ABK37082.1\nMKTLVAILA\n"
    
    # Check if the generated file content matches the expected FASTA output
    expect_equal(fasta_content, expected_fasta)
    
    # Call the function with reduced mode enabled
    result <- convertAlignment2FA(
        aln_file = aln_file,
        lin_file = lin_file,
        fa_outpath = NULL,
        reduced = TRUE
    )
    
    # In reduced mode, we expect one sequence per unique lineage
    expected_fasta <- ">KPhylum_NASpe_ACU53894.1\nMKTIIALSYI\n>KPhylum_NASpe_APJ14606.1\nGSLKVSISY\n>KPhylum_NASpe_ABK37082.1\nMKTLVAILA\n"
    
    # Check if the result matches the expected FASTA output with 
    # reduced sequences
    expect_equal(result, expected_fasta)
    
    # Clean up mock files after test
    file.remove(aln_file)
    file.remove(lin_file)
    
    # Mock data for acc2name
    acc2name <- tibble::tibble(
        AccNum = c("Acc1", "Acc2", "Acc3", "Acc4"),
        Name = c("Name1", "Name2", "Name3", "Name4")
    )
    
    # Sample FASTA line with Acc1
    line <- ">Acc1 AAAATTTTCCCCGGGG"
    
    # Call the function
    result <- mapAcc2Name(line, acc2name)
    
    # Expected output
    expected <- ">Name1"
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Sample FASTA line with an unknown accession number
    line <- ">Acc5 TTTTCCCCGGGGAAAA"
    
    # Call the function
    result <- mapAcc2Name(line, acc2name)
    
    # Expected output (no match found)
    expected <- ">NA"
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Empty line
    line <- ""
    
    # Call the function
    result <- mapAcc2Name(line, acc2name)
    
    # Expected output (should return an empty or NA-like result)
    expect_equal(result, ">NA")
    
    # Sample FASTA line with extra spaces
    line <- ">Acc2    TTTTCCCCGGGGAAAA"
    
    # Call the function
    result <- mapAcc2Name(line, acc2name)
    
    # Expected output
    expected <- ">Name2"
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Sample FASTA line with a lowercase accession number
    # Acc1 is in lowercase, but acc2name has uppercase
    line <- ">acc1 AAAATTTTCCCCGGGG"  
    
    # Call the function
    result <- mapAcc2Name(line, acc2name)
    
    # Expected output: should not find the accession number 
    # because of case difference
    expected <- ">NA"
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Mock data for acc2name
    acc2name <- tibble::tibble(
        AccNum = c("Acc1", "Acc2", "Acc3"),
        Name = c("Name1", "Name2", "Name3")
    )
    
    # Mock FASTA file lines
    mock_fasta <- c(
        ">Acc1 description",
        "AAAATTTTCCCCGGGG",
        ">Acc2 description",
        "GGGGTTTTAAAACCCC",
        ">Acc3 description",
        "CCCCAAAAGGGGTTTT"
    )
    
    # Create a temporary file for the mock FASTA file
    fa_path <- tempfile(fileext = ".fa")
    writeLines(mock_fasta, fa_path)
    
    # Output file path
    outpath <- tempfile(fileext = ".fa")
    
    # Call the renameFA function
    result <- renameFA(fa_path, outpath, 
                       replacement_function = mapAcc2Name, 
                       acc2name = acc2name)
    
    # Expected result after renaming
    expected <- c(
        ">Name1",
        "AAAATTTTCCCCGGGG",
        ">Name2",
        "GGGGTTTTAAAACCCC",
        ">Name3",
        "CCCCAAAAGGGGTTTT"
    )
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Call the renameFA function
    result <- renameFA(fa_path, outpath, 
                       replacement_function = mapAcc2Name, 
                       acc2name = acc2name)
    
    # Extract the sequence lines (should not change)
    sequence_lines <- result[c(2, 4, 6)]
    expected_sequences <- c("AAAATTTTCCCCGGGG", 
                            "GGGGTTTTAAAACCCC", 
                            "CCCCAAAAGGGGTTTT")
    
    # Check if sequence lines match the original content
    expect_equal(sequence_lines, expected_sequences)
    
    # Custom function that just appends '_custom' to the accession number
    custom_replacement <- function(line, ...) {
        acc <- str_sub(line, 2) # Remove '>'
        return(paste0(">", acc, "_custom"))
    }
    
    # Call the renameFA function with the custom function
    result <- renameFA(fa_path, outpath, 
                       replacement_function = custom_replacement)
    
    # Expected result after renaming with the custom function
    expected <- c(
        ">Acc1 description_custom",
        "AAAATTTTCCCCGGGG",
        ">Acc2 description_custom",
        "GGGGTTTTAAAACCCC",
        ">Acc3 description_custom",
        "CCCCAAAAGGGGTTTT"
    )
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Create an empty file
    empty_fa <- tempfile(fileext = ".fa")
    writeLines(character(0), empty_fa)
    
    # Call the renameFA function with the empty file
    result <- renameFA(empty_fa, outpath, 
                       replacement_function = mapAcc2Name, 
                       acc2name = acc2name)
    
    # Expected output is an empty character vector
    expect_equal(result, NULL)
    
    # Call the renameFA function with multiple sequences
    result <- renameFA(fa_path, outpath, 
                       replacement_function = mapAcc2Name, 
                       acc2name = acc2name)
    
    # Expected result
    expected <- c(
        ">Name1",
        "AAAATTTTCCCCGGGG",
        ">Name2",
        "GGGGTTTTAAAACCCC",
        ">Name3",
        "CCCCAAAAGGGGTTTT"
    )
    
    # Check if the result matches the expected output
    expect_equal(result, expected)
    
    # Setup: Create a temporary directory for alignment files
    aln_dir <- tempdir()
    fa_dir <- tempdir()
    lin_dir <- tempdir()
    
    # Create tabular alignment files (two columns: AccNum <TAB> Alignment)
    aln1 <- 
      tibble::tibble(
        AccNum = c("Acc1", "Acc2"),
        Alignment = c("ATGCNNNN", "ATGCNNNN")
        )

    aln2 <- 
      tibble::tibble(
        AccNum = c("Acc3", "Acc4"),
        Alignment = c("ATGCNNNN", "ATGCNNNN")
        )

    # Write without column names so read_tsv(..., col_names = c("AccNum","Alignment")) works
    readr::write_tsv(aln1, file.path(aln_dir, "test1.aln"), col_names = FALSE)
    readr::write_tsv(aln2, file.path(aln_dir, "test2.aln"), col_names = FALSE)
    
    # Create mock lineage data (AccNum, Species, Lineage)
    mock_lineage_data <- tibble::tibble(
        AccNum = c("Acc1", "Acc2", "Acc3", "Acc4"),
        Species = c("Species A", "Species B", "Species C", "Species D"),
        Lineage = c("Kingdom>Phylum1", "Kingdom>Phylum1", 
                    "Kingdom>Phylum2", "Kingdom>Phylum2")
    )
    
    # Save the mock lineage data as a TSV file
    write_tsv(mock_lineage_data, file.path(lin_dir, "all_semiclean.txt"))
    
    # Call the function
    generateAllAlignments2FA(aln_path = aln_dir, 
                             fa_outpath = fa_dir, 
                             lin_file = file.path(lin_dir, "all_semiclean.txt"))
    
    # Check if output files are created
    output_files <- list.files(fa_dir, pattern = "*.fa", full.names = TRUE)
    expect_true(length(output_files) > 0)
    
    # Create a temporary empty directory and test for error handling
    empty_dir <- tempdir()
    fa_dir <- tempdir()

    # Setup: Create a temporary directory for alignment files
    aln_dir <- tempdir()
    fa_dir <- tempdir()
    
  # Create mock alignment data as a tibble
    mock_aln <- tibble(
        AccNum    = c("Acc1", "Acc2"),
        Alignment = c("ATGCNNNN", "ATGCNNNN")
    )

    # Write as a tab-delimited file without header row
    # (so read_tsv(..., col_names = c("AccNum","Alignment")) will work correctly)
    write_tsv(mock_aln,
            file = file.path(aln_dir, "test1.aln"),
            col_names = FALSE)
    # # Create mock alignment files
    # writeLines(c(">Acc1\nATGCNNNN\n>Acc2\nATGCNNNN\n"), 
    #            file.path(aln_dir, "test1.aln"))
    
    # Create mock lineage data
    mock_lineage_data <- tibble::tibble(
        AccNum = c("Acc1", "Acc2"),
        Species = c("Species A", "Species B"),
        Lineage = c("Kingdom>Phylum1", "Kingdom>Phylum1")
    )
    
    # Save the mock lineage data as a TSV file
    write_tsv(mock_lineage_data, file.path(lin_dir, "all_semiclean.txt"))
    
    # Call the function with reduced = TRUE
    generateAllAlignments2FA(aln_path = aln_dir, 
                             fa_outpath = fa_dir, 
                             lin_file = file.path(lin_dir, "all_semiclean.txt"), 
                             reduced = TRUE)
    
    # Check if output files are created
    output_files <- list.files(fa_dir, pattern = "*.fa", full.names = TRUE)
    expect_true(length(output_files) > 0)
    
    # Check the content of the output file
    fasta_content <- readLines(outpath)
    expect_equal(fasta_content[1], ">Name1")
    
    # Sample Data for Testing
    prot_data <- tibble::tibble(
        AccNum = c("A001", "A002", "A003", "A004", "A005"),
        Lineage = c("L1", "L1", "L2", "L2", NA)
    )
    
    result <- createRepresentativeAccNum(prot_data, reduced = "Lineage", accnum_col = "AccNum")
    
    expect_equal(result, c("A001", "A003"))  # A001 for L1 and A003 for L2
    
    result <- createRepresentativeAccNum(prot_data, reduced = "Lineage", accnum_col = "AccNum")
    
    expect_false(any(is.na(result)))  # No NA values in result
    expect_equal(length(result), 2)    # Still only 2 unique Lineages
    
    prot_data_with_duplicates <- tibble::tibble(
        AccNum = c("A001", "A001", "A002", "A003", "A004", "A004"),
        Lineage = c("L1", "L1", "L2", "L2", "L3", "L3")
    )
    
    result <- createRepresentativeAccNum(prot_data_with_duplicates, 
                                    reduced = "Lineage", accnum_col = "AccNum")
    
    # A001 for L1, A002 for L2, and one for L3
    expect_equal(result, c("A001", "A002", "A004"))  
    
    empty_data <- tibble::tibble(AccNum = character(), Lineage = character())
    
    result <- createRepresentativeAccNum(empty_data, 
                                    reduced = "Lineage", 
                                    accnum_col = "AccNum")
    
    expect_equal(result, NULL)  # Expecting an empty character vector
    
    # Create a mock FASTA file for testing
    create_mock_fasta <- function(file_path) {
        seqs <- c(">seq1\nACGT\n>seq2\nACGT\n>seq3\nACGT")
        writeLines(seqs, con = file_path)
    }
    
    # Path for the mock FASTA file
    mock_fasta_path <- "mock_sequences.fasta"
    create_mock_fasta(mock_fasta_path)
    
    aligned <- alignFasta(mock_fasta_path, tool = "Muscle")
    
    # Check the class of the returned object
    expect_s4_class(aligned, "MsaAAMultipleAlignment")
    
    aligned <- alignFasta(mock_fasta_path, tool = "ClustalO")
    
    # Check the class of the returned object
    expect_s4_class(aligned, "MsaAAMultipleAlignment")
    
    aligned <- alignFasta(mock_fasta_path, tool = "ClustalW")
    
    # Check the class of the returned object
    expect_s4_class(aligned, "MsaAAMultipleAlignment")
    
    # Expect error for non-existing file
    expect_error(alignFasta("non_existing_file.fasta", tool = "Muscle"), 
                 "Error: The FASTA file does not exist: ") 
    
    # Clean up mock fasta file after tests
    file.remove(mock_fasta_path)
})
