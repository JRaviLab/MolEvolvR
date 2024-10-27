# config.R
url <- "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"
destfile <- "broadstreet-v3.3.0.tar.bz2"

# Download the file
download.file(url, destfile)

#Extract the file
if (!require("R.utils")) {
  install.packages("R.utils")
  library(R.utils)
}


# Extract the tar file
untar("broadstreet-v3.3.0.tar", exdir = "CARD_data")


# Map CARD Short Name
# Install and Load dplyr and readr
packages <- c("dplyr", "readr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

# Parse the required files using readr::read_delim
aro_index <- read_delim("CARD_data/aro_index.tsv", delim = "\t", col_names = TRUE)
antibiotics_data <- read_delim("CARD_data/shortname_antibiotics.tsv", delim = "\t", col_names = TRUE)
pathogens_data <- read_delim("CARD_data/shortname_pathogens.tsv", delim = "\t", col_names = TRUE)

# Extract pathogen, gene, drug, and include Protein.Accession from 'CARD.Short.Name'
aro_index_clean <- aro_index %>%
  mutate(
    pathogen = sapply(strsplit(CARD.Short.Name, "_"), `[`, 1),  # Extract pathogen
    gene = sapply(strsplit(CARD.Short.Name, "_"), `[`, 2),      # Extract gene
    drug = ifelse(sapply(strsplit(CARD.Short.Name, "_"), length) == 3,   # Extract drug
                  sapply(strsplit(CARD.Short.Name, "_"), `[`, 3), NA),
    Protein.Accession = Protein.Accession  # Include the Protein.Accession column
  )

# Extract pathogen, gene, drug, and include Protein.Accession from 'CARD Short Name'
library(dplyr)
library(purrr)
library(stringr)

# Extract pathogen, gene, drug, and include Protein.Accession from 'CARD Short Name'
resistance_profile <- aro_index %>%
  mutate(
    split_name = strsplit(`CARD Short Name`, "_"),  # Split the CARD Short Name
    pathogen = map_chr(split_name, ~ if (length(.) >= 1) .[1] else NA),  # Extract pathogen
    gene = map_chr(split_name, ~ if (length(.) >= 2) .[2] else NA),      # Extract gene
    drug = map_chr(split_name, ~ if (length(.) == 3) .[3] else NA),      # Extract drug
    Protein.Accession = `Protein Accession`  # Include the Protein Accession column
  ) %>%
  select(-split_name)

# View the cleaned data
head(resistance_profile)


# Merge resistance_profile with the antibiotics_data and pathogens_data
# For merging with antibiotics_data
merged_data_antibiotics <- left_join(resistance_profile, antibiotics_data,
                                     by = c("drug" = "AAC Abbreviation"),
                                     relationship = "many-to-many")

# View the merged data
head(merged_data_antibiotics)


# For merging with pathogens_data
merged_data_pathogens <- left_join(merged_data_antibiotics, pathogens_data,
                                   by = c("pathogen" = "Abbreviation"))

# View the resulting merged data
head(merged_data_pathogens)


#filter out rows where pathogen is empty
cleaned_data <- merged_data_pathogens %>%
  distinct() %>%
  filter(!is.na(Pathogen))
View(cleaned_data)

# Group by Pathogen, Gene, Drug, and Protein_Accession, then summarize Antibiotic information
summarized_data <- cleaned_data %>%
  group_by(Pathogen = Pathogen, Gene = gene, Drug = drug, Protein_Accession = Protein_Accession) %>%
  summarize(Antibiotic_Info = paste(unique(Molecule), collapse = ", ")) %>%
  arrange(Pathogen, Gene, Drug, Protein_Accession)

# Filter for Staphylococcus aureus and DAP (Bug-Drug of Interest)
staph_aureus_dap_combinations <- summarized_data %>%
  filter(Pathogen == "Staphylococcus aureus", Drug == "DAP")

# View the filtered data
head(staph_aureus_dap_combinations)

# FASTA sequences
# Install and Load required packages
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
if (!requireNamespace("XML", quietly = TRUE)) {
  install.packages("XML")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

library(rentrez)
library(XML)
library(stringr)

# Fetch FASTA sequence from Entrez
fetch_fasta_sequence <- function(Protein_Accession) {
  tryCatch({
    # Fetch the FASTA sequence using Entrez
    fasta_seq <- rentrez::entrez_fetch(db = "protein",
                                       id = Protein_Accession,
                                       rettype = "fasta",
                                       retmode = "text")

    if (!is.null(fasta_seq)) {
      # Ensure the first line starts with ">"
      if (!grepl("^>", fasta_seq[1])) {
        fasta_seq[1] <- paste0(">", fasta_seq[1])
      }

      # Split the sequence into lines
      lines <- str_split(fasta_seq, "\n")[[1]]

      # Join the lines back together
      fasta_seq <- paste(lines, collapse = "\n")

      return(fasta_seq)
    } else {
      warning(paste("Failed to retrieve FASTA sequence for protein accession:", Protein_Accession))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Error fetching FASTA sequence for protein accession:", Protein_Accession, ":", e$message))
    return(NULL)
  })
}


# Loop through staph_aureus_dap_combinations to fetch and save FASTA sequences
combined_sequences <- character()

for (i in 1:nrow(staph_aureus_dap_combinations)) {
  # Fetch FASTA sequence for each protein accession
  protein_accession <- staph_aureus_dap_combinations$Protein_Accession[i]
  fasta_sequence <- fetch_fasta_sequence(protein_accession)

  if (!is.null(fasta_sequence)) {
    combined_sequences <- c(combined_sequences, fasta_sequence)
  }
}

# Save the combined FASTA sequences
filename <- "Staph_aureus_Daptomycin_sequences.fasta"

writeLines(combined_sequences, filename)

# Read the FASTA file
fasta_content <- readLines(filename)

# Display the contents
cat(fasta_content, sep = "\n")
