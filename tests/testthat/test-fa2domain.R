context("fa2domain")
test_that("fa2domain", {
    library(mockery)
    library(readr)
    # runIPRScan
    # Define file paths using system.file to locate files in package
    filepath_fasta <- system.file("tests", "example_fasta.fa", 
                                  package = "MolEvolvR")
    filepath_out <- tempfile()  # Temporary file for output
    
    # Set application options
    mock_appl_single <- "Pfam"
    mock_appl_multiple <- c("Pfam", "Gene3D")
    
    # Create a sample TSV file in extdata and read it
    sample_tsv_path <- system.file("tests", "example_iprscan_valid.tsv", 
                                   package = "MolEvolvR")
    # Assumes TSV format
    sample_tsv <- read.csv(sample_tsv_path, sep = "\t", header = TRUE) 
    
    # Mock the system function to avoid running real command
    mock_system <- mock(0L)  # Simulate successful system call
    
    # Patch the system and readIPRScanTSV functions
    stub(runIPRScan, "system", mock_system)
    stub(runIPRScan, "readIPRScanTSV", function(x) read.csv(sample_tsv_path, sep = "\t"))
    
    ## TEST 1: Command construction for single application
    runIPRScan(filepath_fasta, filepath_out, appl = mock_appl_single)
    expected_cmd_single <- glue::glue("iprscan -i {filepath_fasta} -b 
                                      {filepath_out} --cpu 4 -f TSV --appl 
                                      {mock_appl_single}")
    
    # Verify system command was called with expected string
    expect_called(mock_system, 1)
    expect_call(mock_system, 1, system(expected_cmd_single))
    
    ## TEST 2: Command construction for multiple applications
    runIPRScan(filepath_fasta, filepath_out, appl = mock_appl_multiple)
    expected_cmd_multiple <- glue::glue("iprscan -i {filepath_fasta} -b 
                                        {filepath_out} --cpu 4 -f TSV --appl 
                                        {paste(mock_appl_multiple, collapse=',')}")
    
    expect_called(mock_system, 2)
    expect_call(mock_system, 2, system(expected_cmd_multiple))
    
    ## TEST 3: Real result from reading TSV file
    result <- runIPRScan(filepath_fasta, filepath_out, appl = mock_appl_single)
    
    # Check if the result matches the expected data from the TSV file
    expect_equal(result, sample_tsv)
    
    ## TEST 4: Error handling when system command fails
    mock_system_fail <- mock(1L)  # Simulate non-zero exit code
    stub(runIPRScan, "system", mock_system_fail)
    
    # Expect a warning and return NULL on failure
    expect_warning(result_fail <- runIPRScan(filepath_fasta, filepath_out, 
                                             appl = mock_appl_single),
                   "interproscan exited with non-zero code")
    expect_null(result_fail)
    
    ## TEST 5: Error handling for missing or invalid inputs
    # Invalid `filepath_fasta`
    expect_error(runIPRScan(NULL, filepath_out, appl = mock_appl_single), 
                 "filepath_fasta cannot be NULL or empty")
    
    # Invalid `filepath_out`
    expect_error(runIPRScan(filepath_fasta, NULL, appl = mock_appl_single), 
                 "filepath_out cannot be NULL or empty")
    
    # Invalid `appl`
    expect_error(runIPRScan(filepath_fasta, filepath_out, appl = "InvalidApp"), 
                 "Invalid application specified")
    
    # readIPRScanTSV
    # Read the TSV file using the function
    df_ipr <- readIPRScanTSV(sample_tsv_path)
    
    # Check that the returned object is a data frame
    expect_s4_class(df_ipr, "data.frame")
    
    # getIPRScanColNames
    # Call the function to get the column names
    col_names <- getIPRScanColNames()
    
    # Check that the result is a character vector
    expect_type(col_names, "character")
    
    # Define the expected column names
    expected_col_names <- c(
        "AccNum", "SeqMD5Digest", "SLength", "Analysis",
        "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
        "Status", "RunDate", "IPRAcc", "IPRDesc"
    )
    
    # Check that the column names match exactly
    expect_equal(col_names, expected_col_names)
    expect_type(col_names, "character")
    
    # Ensure there are exactly 13 columns
    expect_length(col_names, 13)
    
    # getIPRScanColTypes
    col_types <- getIPRScanColTypes()
    
    # Check that col_types is of the expected class
    # readr::cols() returns col_spec object
    expect_s4_class(col_types, "col_spec")  
    
    # Verify that each column has the correct type
    expect_equal(col_types$cols$AccNum, col_character())
    expect_equal(col_types$cols$SeqMD5Digest, col_character())
    expect_equal(col_types$cols$SLength, col_integer())
    expect_equal(col_types$cols$Analysis, col_character())
    expect_equal(col_types$cols$DB.ID, col_character())
    expect_equal(col_types$cols$SignDesc, col_character())
    expect_equal(col_types$cols$StartLoc, col_integer())
    expect_equal(col_types$cols$StopLoc, col_integer())
    expect_equal(col_types$cols$Score, col_double())
    expect_equal(col_types$cols$Status, col_character())
    expect_equal(col_types$cols$RunDate, col_character())
    expect_equal(col_types$cols$IPRAcc, col_character())
    expect_equal(col_types$cols$IPRDesc, col_character())
    
    # Optionally, check that there are no additional columns defined
    expect_length(col_types$cols, 13)
    
    # createIPRScanDomainTable
    
    # Load the sample FASTA file
    fasta <- Biostrings::readAAStringSet(filepath_fasta)
    
    # Read the sample InterProScan TSV file
    df_iprscan <- readIPRScanTSV(sample_tsv_path)
    
    # Example accession number for testing
    accnum <- df_iprscan$AccNum[1]
    
    # Test case 1: Valid inputs
    df_iprscan_domains <- createIPRScanDomainTable(accnum, fasta, df_iprscan)
    
    # Check that the output is a data frame
    expect_s4_class(df_iprscan_domains, "data.frame")
    
    # Validate the structure of the output
    expect_true(all(c("AccNum", "DB.ID", "StartLoc", "StopLoc", "seq_domain", 
                      "id_domain") %in% names(df_iprscan_domains)))
    
    # Validate the content of the seq_domain column
    # Ensure no empty sequences
    expect_true(all(nchar(df_iprscan_domains$seq_domain) > 0))  
    
    # Validate the id_domain structure
    expect_true(all(grepl("^\\w+-\\w+-\\d+_\\d+$", df_iprscan_domains$id_domain)))
    
    # Test case 2: No matching accession number
    empty_df <- createIPRScanDomainTable("non_existent_accnum", fasta, df_iprscan)
    expect_s3_class(empty_df, "data.frame")
    expect_equal(nrow(empty_df), 0)
    
    # Test case 3: No domains in input data frame
    empty_iprscan <- df_iprscan[0, ]  # Create an empty df_iprscan
    empty_domains_df <- createIPRScanDomainTable(accnum, fasta, empty_iprscan)
    expect_s3_class(empty_domains_df, "data.frame")
    expect_equal(nrow(empty_domains_df), 0)
    
    # convertIPRScanDomainTable2FA

    # Test case 1: Valid domain data
    fasta_domains <- convertIPRScanDomainTable2FA(df_iprscan_domains)
    
    # Check that the output is an AAStringSet
    expect_s4_class(fasta_domains, "AAStringSet")
    
    # Check that the correct number of sequences are returned
    expect_equal(length(fasta_domains), nrow(df_iprscan_domains))
    
    # Check that the names of the sequences match the id_domain column
    expect_equal(names(fasta_domains), df_iprscan_domains$id_domain)
    
    # Test case 2: Empty input data frame
    empty_domains <- convertIPRScanDomainTable2FA(data.frame())
    expect_s4_class(empty_domains, "AAStringSet")
    expect_equal(length(empty_domains), 0)
    
    # Test case 3: Data frame with no domains
    empty_df_iprscan <- df_iprscan[0, ]  # Create an empty df_iprscan
    empty_domains_df <- convertIPRScanDomainTable2FA(empty_df_iprscan)
    expect_s4_class(empty_domains_df, "AAStringSet")
    expect_equal(length(empty_domains_df), 0)
    
    # getDomainsFromFA
    # Test case 1: Valid input
    fasta_domains <- getDomainsFromFA(fasta, df_iprscan)
    
    # Check that the output is an AAStringSet
    expect_s4_class(fasta_domains, "AAStringSet")
    
    # Check that the output contains the expected sequences
    expect_true(length(fasta_domains) > 0)  # Ensure there are some domains extracted
    
    # Test case 2: Empty input FASTA
    empty_fasta <- Biostrings::AAStringSet()
    empty_fasta_domains <- getDomainsFromFA(empty_fasta, df_iprscan)

    expect_s4_class(empty_fasta_domains, "AAStringSet")
    expect_equal(length(empty_fasta_domains), 0)
    
    # Test case 3: Empty input df_iprscan
    empty_iprscan <- data.frame()  # Create an empty df_iprscan
    empty_domains_iprscan <- getDomainsFromFA(fasta, empty_iprscan)
    
    expect_s4_class(empty_domains_iprscan, "AAStringSet")
    expect_equal(length(empty_domains_iprscan), 0)
    
    # Test case 4: Verbose output
    expect_warning(
        getDomainsFromFA(fasta, empty_iprscan, verbose = TRUE),
        regexp = "had no domains for the selected analyses"
    )
    
    # Test case 5: Verbose output for some valid accession numbers
    fasta_domains_verbose <- getDomainsFromFA(fasta, df_iprscan, verbose = TRUE)
    
    # Check that the output is still an AAStringSet
    expect_s4_class(fasta_domains_verbose, "AAStringSet")
    
})