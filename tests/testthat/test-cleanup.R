context("cleanup")
test_that("cleanup", {
    # cleanup
    # Test with normal string
    expect_equal(cleanString("Hello World"), "Hello_World")
    
    # Test with multiple spaces
    expect_equal(cleanString("Hello     World"), "Hello_World")
    
    # Test with special characters
    expect_equal(cleanString("Hello@World!"), "HelloWorld")
    
    # Test with alphanumeric characters and underscore
    expect_equal(cleanString("Test_String 123"), "Test_String_123")
    
    # Test with dots
    expect_equal(cleanString("Version 1.0.0"), "Version_1.0.0")
    
    # Test with empty string
    expect_equal(cleanString(""), "")
    
    # Test with only spaces
    expect_equal(cleanString("    "), "_")
    
    # Test with non-alphanumeric characters
    expect_equal(cleanString("~!@#$%^&*()"), "")
    
    # Test with mixed characters
    expect_equal(cleanString("Hello !@#$% World"), "Hello__World")
    
    # Test with trailing and leading spaces
    expect_equal(cleanString("  Test  "), "_Test_")
    
    # Test with numbers and underscores
    expect_equal(cleanString("Name_123 Test"), "Name_123_Test")
    
    # extractAccNum
    # Test with a string containing a pipe character
    expect_equal(extractAccNum("ID|ABC1234 Some Description"), "ABC1234")
    
    # Test with a string containing multiple spaces
    expect_equal(extractAccNum("ID|DEF5678    More Info"), "DEF5678")
    
    # Test with a string without a pipe character
    expect_equal(extractAccNum("ABC9876 Some Description"), "ABC9876")
    
    # Test with a string that has leading spaces
    expect_equal(extractAccNum("   ID|GHI1357 Description"), "GHI1357")
    
    # Test with a string that has trailing spaces
    expect_equal(extractAccNum("ID|JKL2468 Description   "), "JKL2468")
    
    # Test with only an accession number
    expect_equal(extractAccNum("XYZ1234"), "XYZ1234")
    
    # Test with a string with only spaces
    expect_equal(extractAccNum("    "), "")
    
    # Test with a string that contains special characters
    expect_equal(extractAccNum("ID|MNO5678_Extra Info"), "MNO5678_Extra")
    
    # ensureUniqAccNum
    # Test with unique accession numbers
    accnums1 <- c("ABC1234", "DEF5678", "GHI9012")
    expect_equal(ensureUniqAccNum(accnums1), c("ABC1234_1", "DEF5678_1", 
                                               "GHI9012_1"))
    
    # Test with duplicate accession numbers
    accnums2 <- c("ABC1234", "ABC1234", "DEF5678", "DEF5678", "GHI9012")
    expect_equal(ensureUniqAccNum(accnums2), 
                 c("ABC1234_1", "ABC1234_2", "DEF5678_1", 
                   "DEF5678_2", "GHI9012_1"))
    
    # Test with all identical accession numbers
    accnums3 <- c("XYZ9999", "XYZ9999", "XYZ9999")
    expect_equal(ensureUniqAccNum(accnums3), 
                 c("XYZ9999_1", "XYZ9999_2", "XYZ9999_3"))
    
    # Test with empty input
    accnums4 <- character(0)
    expect_equal(ensureUniqAccNum(accnums4), character(0))
    
    # Test with a single accession number
    accnums5 <- c("SINGLE_ACC")
    expect_equal(ensureUniqAccNum(accnums5), c("SINGLE_ACC_1"))
    
    # Test with mixed duplicate and unique accession numbers
    accnums6 <- c("A", "B", "A", "C", "B", "B")
    expect_equal(ensureUniqAccNum(accnums6), 
                 c("A_1", "A_2", "B_1", "B_2", "B_3", "C_1"))
    
    # cleanFAHeaders
    fasta_sample <- c(
        ">sp|P12345|ProteinA Description 1",
        ">sp|P67890|ProteinB Description 2",
        ">sp|P12345|ProteinA Description 3",
        ">sp|P67890|ProteinB Description 4"
    )
    names(fasta_sample) <- fasta_sample  # Set names to headers
    
    # Run the function
    cleaned_fasta <- cleanFAHeaders(fasta_sample)
    
    # Expected headers after processing
    expected_headers <- c("P12345_1", "P12345_2", "P67890_1", "P67890_2")
    
    # Check if the names of cleaned_fasta match expected_headers
    expect_equal(names(cleaned_fasta), expected_headers)
    
    # Check that the contents of cleaned_fasta remain unchanged
    expect_equal(as.vector(cleaned_fasta), as.vector(fasta_sample))
    
    fasta_unique <- c(
        ">sp|P12345|UniqueProteinA",
        ">sp|P67890|UniqueProteinB"
    )
    names(fasta_unique) <- fasta_unique
    
    cleaned_unique_fasta <- cleanFAHeaders(fasta_unique)
    
    expected_unique_headers <- c("P12345_1", "P67890_1")
    expect_equal(names(cleaned_unique_fasta), expected_unique_headers)
    
    # Sample input data
    prot_data <- tibble::tibble(
        DomArch = c("ABC123", "-", "NA", "", "XYZ789", " "),
        other_col = c(1, 2, 3, 4, 5, 6)
    )
    
    # Expected output after removing rows
    expected_output <- tibble::tibble(
        DomArch = c("ABC123", "XYZ789"),
        other_col = c(1, 5)
    )
    
    # Run the function
    result <- removeEmptyRows(prot_data, by_column = "DomArch")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Case 1: No rows removed
    prot_data_no_removal <- tibble::tibble(
        DomArch = c("ABC123", "XYZ789"),
        other_col = c(1, 2)
    )
    expect_equal(removeEmptyRows(prot_data_no_removal), prot_data_no_removal)
    
    # Case 2: All rows removed
    prot_data_all_removed <- tibble::tibble(
        DomArch = c("-", "NA", "", " "),
        other_col = c(1, 2, 3, 4)
    )
    expect_equal(removeEmptyRows(prot_data_all_removed), 
                 tibble::tibble(DomArch = character(0), other_col = numeric(0)))
    
    # Case 3: Empty data frame
    prot_data_empty <- tibble::tibble(DomArch = character(0), 
                                      other_col = numeric(0))
    expect_equal(removeEmptyRows(prot_data_empty), prot_data_empty)
    
    # Input data with repeated domains
    prot_data <- tibble::tibble(
        DomArch = c("A B B C", "X X Y", "P P P Q", "R R R S"),
        other_col = c(1, 2, 3, 4)
    )
    
    # Input data with repeated domains
    prot_data <- tibble::tibble(
        DomArch = c("A A A", "B B", "C C C D D"),
        other_col = c(1, 2, 3)
    )
    
    # Input data with repeated and single question marks
    prot_data <- tibble::tibble(
        GenContext = c("A ? ? B", "? ?", "C ?? C", "D ? > ? D"),
        other_col = c(1, 2, 3, 4)
    )
    
    # Input data with single question marks only
    prot_data <- tibble::tibble(
        GenContext = c("?", "? ? ?", "A ? B"),
        other_col = c(1, 2, 3)
    )
    
    # Expected output after replacing single question marks
    expected_output <- tibble::tibble(
        GenContext = c("X", "X(s)", "A X B"),
        other_col = c(1, 2, 3)
    )
    
    # Run the function
    result <- replaceQuestionMarks(prot_data, by_column = "GenContext")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Input data containing asterisks
    query_data <- tibble::tibble(
        GenContext = c("A * B", "*C*D*", "E*F*"),
        other_col = c(1, 2, 3)
    )
    
    # Expected output after removing asterisks
    expected_output <- tibble::tibble(
        GenContext = c("A  B", "CD", "EF"),
        other_col = c(1, 2, 3)
    )
    
    # Run the function
    result <- removeAsterisks(query_data, colname = "GenContext")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Input data with no asterisks
    query_data <- tibble::tibble(
        GenContext = c("A B", "C D", "E F"),
        other_col = c(1, 2, 3)
    )
    
    # Expected output (no changes)
    expected_output <- query_data
    
    # Run the function
    result <- removeAsterisks(query_data, colname = "GenContext")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data
    prot <- tibble::tibble(
        DomArch = c("A", "B", "A", "C", "D"),
        value = c(1, 2, 3, 4, 5)
    )
    
    # Expected output after removing rows where `DomArch` appears only once
    expected_output <- tibble::tibble(
        DomArch = c("A", "A"),
        value = c(1, 3)
    )
    
    # Run the function
    result <- removeTails(prot, by_column = "DomArch", keep_domains = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Input data with no single occurrence rows
    prot <- tibble::tibble(
        DomArch = c("A", "A", "B", "B"),
        value = c(1, 2, 3, 4)
    )
    
    # Expected output (should remain unchanged)
    expected_output <- prot
    
    # Run the function
    result <- removeTails(prot, by_column = "DomArch", keep_domains = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with special characters and extra spaces in the species names
    prot <- tibble::tibble(
        Species.orig = c("Escherichia coli sp.", 
                         "Bacillus str. subtilis", 
                         "Lactobacillus = plantarum", 
                         "Staphylococcus aureus"),
        value = c(1, 2, 3, 4)
    )
    
    # Expected output after cleaning species names
    expected_output <- tibble::tibble(
        Species.orig = c("Escherichia coli sp.", 
                         "Bacillus str. subtilis", 
                         "Lactobacillus = plantarum", 
                         "Staphylococcus aureus"),
        value = c(1, 2, 3, 4),
        Species = c("Escherichia coli sp", 
                    "Bacillus str subtilis", 
                    "Lactobacillus plantarum", 
                    "Staphylococcus aureus")
    )
    
    # Run the function
    result <- cleanSpecies(prot, removeEmptyRows = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with an empty Species entry
    prot <- tibble::tibble(
        Species.orig = c("Escherichia coli sp.", "", 
                         "Lactobacillus = plantarum", 
                         "Staphylococcus aureus"),
        value = c(1, 2, 3, 4)
    )
    
    # Expected output after cleaning and removing empty rows
    expected_output <- tibble::tibble(
        Species.orig = c("Escherichia coli sp.", 
                         "Lactobacillus = plantarum", 
                         "Staphylococcus aureus"),
        value = c(1, 3, 4),
        Species = c("Escherichia coli sp", 
                    "Lactobacillus plantarum", 
                    "Staphylococcus aureus")
    )
    
    # Run the function with removeEmptyRows = TRUE
    result <- cleanSpecies(prot, removeEmptyRows = TRUE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with original ClustName
    prot <- tibble::tibble(
        ClustName.orig = c("SIG+TM+TM", "ABC+DEF", "XYZ+SIG", "TM+TM"),
        value = c(1, 2, 3, 4)
    )
    
    # Domains to rename
    domains_rename <- tibble::tibble(
        old = c("SIG", "ABC"),
        new = c("Signal", "ABC_Transporter")
    )
    
    # Domains to keep
    domains_keep <- tibble::tibble(
        domains = c("Signal", "ABC_Transporter")
    )
    
    # Expected output after renaming and filtering
    expected_output <- tibble::tibble(
        ClustName.orig = c("SIG+TM+TM", "ABC+DEF", "XYZ+SIG"),
        value = c(1, 2, 3),
        ClustName = c("Signal+TM+TM", "ABC_Transporter+DEF", "XYZ+Signal")
    )
    
    # Run the function
    result <- cleanClusters(prot, domains_rename, domains_keep, 
                            condenseRepeatedDomains = FALSE, 
                            removeTails = FALSE, 
                            removeEmptyRows = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with ClustName containing tails
    prot <- tibble::tibble(
        ClustName.orig = c("SIG+TM+1", "ABC+DEF", "XYZ+SIG+2"),
        value = c(1, 2, 3)
    )
    
    # Domains to rename (empty for this test)
    domains_rename <- tibble::tibble(
        old = character(0),
        new = character(0)
    )
    
    # Domains to keep (empty for this test)
    domains_keep <- tibble::tibble(
        domains = character(0)
    )
    
    # Expected output after removing tails
    expected_output <- tibble::tibble(
        ClustName.orig = c("ABC+DEF"),
        value = c(2),
        ClustName = c("ABC+DEF")
    )
    
    # Sample input data
    prot <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM", "ABC+DEF", "XYZ+SIG", "TM+TM"),
        value = c(1, 2, 3, 4)
    )
    
    # Domains to rename
    domains_rename <- tibble::tibble(
        old = c("SIG", "ABC"),
        new = c("Signal", "ABC_Transporter")
    )
    
    # Domains to keep
    domains_keep <- tibble::tibble(
        domains = c("Signal", "ABC_Transporter")
    )
    
    # Expected output after renaming and filtering
    expected_output <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM", "ABC+DEF", "XYZ+SIG"),
        value = c(1, 2, 3),
        DomArch = c("Signal+TM+TM", "ABC_Transporter+DEF", "XYZ+Signal")
    )
    
    # Run the function
    result <- cleanDomainArchitecture(prot, old = "DomArch.orig", new = "DomArch",
                                      domains_keep = domains_keep,
                                      domains_rename = domains_rename,
                                      condenseRepeatedDomains = FALSE,
                                      removeTails = FALSE,
                                      removeEmptyRows = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with repeated domains
    prot <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM+TM", "ABC+ABC+DEF", "XYZ+SIG+SIG"),
        value = c(1, 2, 3)
    )
    
    # Domains to rename (empty for this test)
    domains_rename <- tibble::tibble(
        old = character(0),
        new = character(0)
    )
    
    # Domains to keep
    domains_keep <- tibble::tibble(
        domains = c("SIG", "ABC")
    )
    
    # Expected output after condensing repeated domains
    expected_output <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM+TM", "ABC+ABC+DEF", "XYZ+SIG+SIG"),
        value = c(1, 2, 3),
        DomArch = c("SIG+TM(s)", "ABC(s)+DEF", "XYZ+SIG(s)")
    )
    
    # Run the function with condenseRepeatedDomains = TRUE
    result <- cleanDomainArchitecture(prot, old = "DomArch.orig", new = "DomArch",
                                      domains_keep = domains_keep,
                                      domains_rename = domains_rename,
                                      condenseRepeatedDomains = TRUE,
                                      removeTails = FALSE,
                                      removeEmptyRows = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with an empty DomArch entry
    prot <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM", "", "ABC+DEF"),
        value = c(1, 2, 3)
    )
    
    # Domains to rename (empty for this test)
    domains_rename <- tibble::tibble(
        old = character(0),
        new = character(0)
    )
    
    # Domains to keep
    domains_keep <- tibble::tibble(
        domains = c("SIG", "ABC")
    )
    
    # Expected output after removing empty rows
    expected_output <- tibble::tibble(
        DomArch.orig = c("SIG+TM+TM", "ABC+DEF"),
        value = c(1, 3),
        DomArch = c("SIG+TM+TM", "ABC+DEF")
    )
    
    # Run the function with removeEmptyRows = TRUE
    result <- cleanDomainArchitecture(prot, old = "DomArch.orig", new = "DomArch",
                                      domains_keep = domains_keep,
                                      domains_rename = domains_rename,
                                      condenseRepeatedDomains = FALSE,
                                      removeTails = FALSE,
                                      removeEmptyRows = TRUE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data with question marks
    prot <- tibble::tibble(
        DomArch.orig = c("SIG+TM???", "ABC+???DEF", "XYZ+SIG"),
        value = c(1, 2, 3)
    )
    
    # Domains to rename (empty for this test)
    domains_rename <- tibble::tibble(
        old = character(0),
        new = character(0)
    )
    
    # Domains to keep
    domains_keep <- tibble::tibble(
        domains = c("SIG", "ABC")
    )
    
    # Expected output after replacing question marks
    expected_output <- tibble::tibble(
        DomArch.orig = c("SIG+TM???", "ABC+???DEF", "XYZ+SIG"),
        value = c(1, 2, 3),
        DomArch = c("SIG+TMXXX", "ABC+XXXDEF", "XYZ+SIG")
    )
    
    # Run the function with question mark replacement
    result <- cleanDomainArchitecture(prot, old = "DomArch.orig", new = "DomArch",
                                      domains_keep = domains_keep,
                                      domains_rename = domains_rename,
                                      condenseRepeatedDomains = FALSE,
                                      removeTails = FALSE,
                                      removeEmptyRows = FALSE)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data
    prot <- tibble::tibble(
        GeneDescription = c("Gene A.", 
                            "Protein B%2C protein C.", 
                            "Enzyme%2C catalytic."),
        value = c(1, 2, 3)
    )
    
    # Expected output after cleaning
    expected_output <- tibble::tibble(
        GeneDescription = c("Gene A.", 
                            "Protein B%2C protein C.", 
                            "Enzyme%2C catalytic."),
        value = c(1, 2, 3),
        GeneDesc = c("Gene A.", 
                     "Protein B, protein C.", 
                     "Enzyme, catalytic.")
    )
    
    # Run the function
    result <- cleanGeneDescription(prot, "GeneDescription")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data
    # mock data for selectLongestDuplicate()
    prot <- tibble::tibble(
        AccNum = c("P001", "P001", "P002", "P003", "P003", "P003", "P004"),
        Description = c(
            "Short peptide",
            "Much longer peptide description text",
            "Single entry",
            "Tiny desc",
            "Medium length description",
            "This one is definitely the longest description for P003",
            "Unique protein"
        ),
        Length = c(120, 150, 80, 50, 70, 95, 200)
        )
    
    # Expected output after selecting longest duplicates
    expected_output <- tibble::tibble(
        AccNum = c("P001", "P002", "P003", "P004"),
        Description = c(
            "Much longer peptide description text",
            "Single entry",
            "This one is definitely the longest description for P003",
            "Unique protein"
        ),
        Length = c(150, 80, 95, 200)
        )
    
    # Run the function
    result <- selectLongestDuplicate(prot, "Description")
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
    # Sample input data
    prot <- tibble::tibble(
        Lineage = c("Bacteria; Firmicutes; Bacilli; Lactobacillales",
                    "Bacteria; Proteobacteria; Gammaproteobacteria",
                    "Archaea; Euryarchaeota; Methanobacteria")
    )
    
    # Rename mapping
    lins_rename <- tibble::tibble(
        old = c("Bacteria", "Firmicutes", "Archaea"),
        new = c("Bacterium", "Firmicute", "Archaean")
    )
    
    # Expected output after renaming
    expected_output <- tibble::tibble(
        Lineage = c("Bacterium; Firmicute; Bacilli; Lactobacillales",
                    "Bacterium; ProteoBacterium; GammaproteoBacterium",
                    "Archaean; Euryarchaeota; MethanoBacterium")
    )
    
    # Run the function
    result <- cleanLineage(prot, lins_rename)
    
    # Check if the result matches the expected output
    expect_equal(result, expected_output)
    
})
