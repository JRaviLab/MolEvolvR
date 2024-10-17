# config.R

# Specify the URL to download the zipped file directly
data_url <- "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"

# Specify the destination directory where the file will be saved
data_path <- "C:/Users/Administrator/Desktop/Process CARD Data/"
bzip_file <- file.path(data_path, "broadstreet-v3.3.0.tar.bz2")
tar_file <- file.path(data_path, "broadstreet-v3.3.0.tar")

# Function to download the file if it doesn't exist
download_file <- function(url, destfile) {
  if (!file.exists(destfile)) {
    download.file(url, destfile, mode = "wb")
    cat("Downloaded:", destfile, "\n")
  } else {
    cat("File already exists:", destfile, "\n")
  }
}

# Call the function to download the file
download_file(data_url, bzip_file)

# Load the R.utils package for decompression
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}
library(R.utils)

# Decompress the file
if (file.exists(bzip_file)) {
  if (!file.exists(tar_file)) {
    bunzip2(bzip_file, destname = tar_file, remove = FALSE)
    cat("Decompressed:", tar_file, "\n")
  } else {
    cat("Tar file already exists:", tar_file, "\n")
  }
} else {
  cat("Error: Zipped file not found:", bzip_file, "\n")
}


# List the contents of the extraction directory
list.files("CARD_data")

# Parse the ARO_index.tsv file using read.delim
aro_index <- read.delim("CARD_data/ARO_index.tsv", header = TRUE, sep = "\t")


# Map CARD Short Name
# Load necessary library
library(dplyr)

#  Read the files
aro_index <- read.delim("CARD_data/aro_index.tsv", sep = "\t", header = TRUE)
antibiotics_data <- read.delim("CARD_data/shortname_antibiotics.tsv", sep = "\t", header = TRUE)
pathogens_data <- read.delim("CARD_data/shortname_pathogens.tsv", sep = "\t", header = TRUE)


# Mutate data
aro_index <- aro_index %>%
  mutate(
    pathogen = sapply(strsplit(CARD.Short.Name, "_"), `[`, 1),   # First part: Pathogen
    gene = sapply(strsplit(CARD.Short.Name, "_"), `[`, 2),       # Second part: Gene
    drug = ifelse(sapply(strsplit(CARD.Short.Name, "_"), length) == 3,   # Third part (if present): Drug
                  sapply(strsplit(CARD.Short.Name, "_"), `[`, 3), NA),
    Protein.Accession = Protein.Accession   # Include existing Protein.Accession column if present
  )

# View the mutated data
head(aro_index)


# Extract pathogen, gene, drug, and include Protein.Accession from 'CARD.Short.Name'
aro_index_clean <- aro_index %>%
  mutate(
    pathogen = sapply(strsplit(CARD.Short.Name, "_"), `[`, 1),  # Extract pathogen
    gene = sapply(strsplit(CARD.Short.Name, "_"), `[`, 2),      # Extract gene
    drug = ifelse(sapply(strsplit(CARD.Short.Name, "_"), length) == 3,   # Extract drug
                  sapply(strsplit(CARD.Short.Name, "_"), `[`, 3), NA),
    Protein.Accession = Protein.Accession  # Include the Protein.Accession column
  )

# Merge aro_index_clean with the antibiotics_data and pathogens_data
# For merging with antibiotics_data
merged_data_antibiotics <- left_join(aro_index_clean, antibiotics_data,
                                     by = c("gene" = "AAC.Abbreviation"))

# For merging with pathogens_data
merged_data_pathogens <- left_join(merged_data_antibiotics, pathogens_data,
                                   by = c("pathogen" = "Abbreviation"))

# View the resulting merged data
head(merged_data_pathogens)


#filter out rows where pathogen is empty
cleaned_data <- merged_data_pathogens %>%
  distinct() %>%
  filter(!is.na(pathogen))

# Group by Pathogen, Gene, Drug, and Protein.Accession, then summarize Antibiotic information
summarized_data <- cleaned_data %>%
  group_by(Pathogen = Pathogen, Gene = gene, Drug = drug, Protein_Accession = Protein.Accession) %>%
  summarize(Antibiotic_Info = paste(unique(Molecule), collapse = ", ")) %>%
  arrange(Pathogen, Gene, Drug, Protein_Accession)

# Filter for Staphylococcus aureus and DAP (Bug-Drug of Interest)
staph_aureus_dap_combinations <- summarized_data %>%
  filter(Pathogen == "Staphylococcus aureus", Drug == "DAP")

# View the filtered data
head(staph_aureus_dap_combinations)


#Fetch FASTA sequences from Entrez using protein accession
#Load required packages
library(rentrez)
library(XML)
library(stringr)


# Fetch FASTA sequence from Entrez
fetch_fasta_sequence <- function(protein_accession) {
  tryCatch({
    # Fetch the FASTA sequence using Entrez
    fasta_seq <- rentrez::entrez_fetch(db = "protein",
                                       id = protein_accession,
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
      warning(paste("Failed to retrieve FASTA sequence for protein accession:", protein_accession))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Error fetching FASTA sequence for protein accession:", protein_accession, ":", e$message))
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
filename <- "Staphylococcus_aureus_DAP_sequences4.fasta"

writeLines(combined_sequences, filename)
print(paste("Combined FASTA sequences saved as", filename))

# Print the first few lines of the combined file
cat(readLines(filename)[1:10], sep = "\n")


#Run FASTA Sequences Through MolEvolvR Web App

library("devtools")

devtools::check()
