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
untar("broadstreet-v3.3.0.tar.bz2", exdir = "CARD_data")


# Map CARD Short Name

# Parse the required files using readr::read_delim
aro_index <- read_delim("CARD_data/aro_index.tsv", delim = "\t", col_names = TRUE)
antibiotics_data <- read_delim("CARD_data/shortname_antibiotics.tsv", delim = "\t", col_names = TRUE)
pathogens_data <- read_delim("CARD_data/shortname_pathogens.tsv", delim = "\t", col_names = TRUE)



# Extract pathogen, gene, drug, and include Protein.Accession from 'CARD Short Name'
extract_card_info <- function(card_short_name, drug_class, `Protein Accession`, `DNA Accession`) {
  # Split the CARD Short Name by underscores
  split_names <- unlist(strsplit(card_short_name, "_"))
  
  # Initialize variables with defaults
  pathogen <- NA
  gene <- NA
  drug <- drug_class  # Default to Drug Class column
  
  # Determine the information based on the split names and patterns
  if (length(split_names) == 1) {
    # Gene only (single part entry)
    gene <- split_names[1]
    pathogen <- "MULTI"  # Assign MULTI as default for pathogen
  } else if (all(toupper(split_names) == split_names)) {
    # Gene complex (all uppercase entries)
    gene <- card_short_name  # Entire entry as gene
    pathogen <- "MULTI"
  } else if (length(split_names) == 2) {
    # Pathogen-Gene scenario
    pathogen <- split_names[1]
    gene <- split_names[2]
  } else if (length(split_names) == 3) {
    # Pathogen-Gene-Drug scenario
    pathogen <- split_names[1]
    gene <- split_names[2]
    drug <- split_names[3]  # Assign drug from the split entry
  }
  
  # If both pathogen and gene are NA, classify as complex gene
  if (is.na(pathogen) && is.na(gene)) {
    gene <- card_short_name  # Assign entire CARD Short Name as gene
    pathogen <- "MULTI"      # Default to MULTI for pathogen
  }
  
  # Handle Protein Accession
  if (is.na(`Protein Accession`) || `Protein Accession` == "") {
    `Protein Accession` <- `DNA Accession`  # Use DNA Accession if Protein Accession is NA
  }
  
  return(list(Pathogen = pathogen, Gene = gene, Drug = drug, Protein_Accession = `Protein Accession`))
}

# Apply the function to the data frame
resistance_profile_data <- aro_index %>%
  mutate(extracted_info = pmap(list(`CARD Short Name`, `Drug Class`, `Protein Accession`, `DNA Accession`),
                               extract_card_info)) %>%
  unnest_wider(extracted_info)

# View the resulting data frame
print(resistance_profile_data)

# Define a relative path for saving the data
output_path <- file.path("CARD_data", "resistance_profile_data.tsv")

# Save resistance_profile_data to the specified path
write_delim(resistance_profile_data, output_path, delim = "\t")

# Load data
resistance_profile_data <- read_delim(output_path, delim = "\t", col_names = TRUE)
antibiotics_data <- read_delim("CARD_data/shortname_antibiotics.tsv", delim = "\t", col_names = TRUE)
pathogens_data <- read_delim("CARD_data/shortname_pathogens.tsv", delim = "\t", col_names = TRUE)


# Merge the extracted resistance profile data with antibiotics_data on Drug
merged_data_antibiotics <- left_join(
  resistance_profile_data,
  antibiotics_data,
  by = c("Drug" = "AAC Abbreviation"), # Adjusting for abbreviations between datasets
  relationship = "many-to-many"
)

# Merge the result with pathogens_data on Pathogen, renaming Pathogen.y to Pathogen_Full_Name
merged_data_pathogens <- left_join(
  merged_data_antibiotics,
  pathogens_data,
  by = c("Pathogen" = "Abbreviation")
) %>%
  rename(Pathogen_Full_Name = Pathogen.y)

# Assign "Multi-species" to Pathogen_Full_Name where Pathogen values are "MULTI"
merged_data_pathogens <- merged_data_pathogens %>%
  mutate(Pathogen_Full_Name = if_else(Pathogen == "MULTI", "Multi-species", Pathogen_Full_Name))


# Assign "Multi-class" to Molecule where Drug values are full names (not abbreviations)
merged_data_pathogens <- merged_data_pathogens %>%
  mutate(Molecule = if_else(grepl(" ", Drug) | grepl("-", Drug), "Multi-class", Molecule))


#FASTA sequences
#Install and Load required packages
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

# Filter for the target drug (DAP) and pathogen (Staphylococcus aureus)
filtered_data <- merged_data_pathogens %>%
  filter(Drug == "DAP", Pathogen_Full_Name == "Staphylococcus aureus")


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


# Define the output file for the FASTA sequences
output_fasta_file <- "Staph_aureus_Daptomycin_sequences.fasta"

# Initialize an empty character vector to store the sequences
combined_sequences <- character()

# Loop through each Protein Accession in the filtered data to fetch sequences
for (i in 1:nrow(filtered_data)) {
  # Get the Protein Accession ID
  protein_accession <- filtered_data$Protein_Accession[i]
  
  cat("Fetching sequence for Protein Accession:", protein_accession, "\n")  # Debugging message
  
  # Fetch the FASTA sequence
  fasta_sequence <- fetch_fasta_sequence(protein_accession)
  
  # If the sequence was fetched successfully, add it to the combined_sequences vector
  if (!is.null(fasta_sequence)) {
    combined_sequences <- c(combined_sequences, fasta_sequence)
    cat("Successfully fetched sequence for:", protein_accession, "\n")
  } else {
    cat("Failed to fetch sequence for:", protein_accession, "\n")
  }
}

# Check if there are any fetched sequences
if (length(combined_sequences) > 0) {
  # Save all fetched sequences to a FASTA file
  writeLines(combined_sequences, output_fasta_file)
  cat("Sequences saved to", output_fasta_file, "\n")
} else {
  cat("No sequences were fetched, so no FASTA file was created.\n")
}

# Read the contents of the file
fasta_contents <- readLines(output_fasta_file)
  
# Print the contents
cat(fasta_contents, sep = "\n")











