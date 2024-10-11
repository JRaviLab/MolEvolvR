# Specify the URL and the destination path where the file will be saved
url <- "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"
destfile <- "broadstreet-v3.3.0.tar.bz2"

# Download the file
download.file(url, destfile)

#Extract the file
install.packages("R.utils")
library(R.utils)

# Decompress the file
bunzip2("broadstreet-v3.3.0.tar.bz2", destname = "broadstreet-v3.3.0.tar", remove = FALSE)

# Extract the tar file
untar("broadstreet-v3.3.0.tar", exdir = "CARD_data")

# List the contents of the extraction directory
list.files("CARD_data")

# Parse the ARO_index.tsv file using read.delim
aro_index <- read.delim("CARD_data/ARO_index.tsv", header = TRUE, sep = "\t")

# View the first few rows of the data
head(aro_index)


# Map CARD Short Name
# Load necessary library
library(dplyr)

#  Read the files
aro_index <- read.delim("CARD_data/aro_index.tsv", sep = "\t", header = TRUE)
antibiotics_data <- read.delim("CARD_data/shortname_antibiotics.tsv", sep = "\t", header = TRUE)
pathogens_data <- read.delim("CARD_data/shortname_pathogens.tsv", sep = "\t", header = TRUE)

# View the first few rows to ensure data is loaded correctly
head(aro_index)
head(antibiotics_data)
head(pathogens_data)

# Split CARD.Short.Name into pathogen, gene, and drug
aro_index <- aro_index %>%
  mutate(
    pathogen = sapply(strsplit(CARD.Short.Name, "_"), `[`, 1),   # First part: Pathogen
    gene = sapply(strsplit(CARD.Short.Name, "_"), `[`, 2),       # Second part: Gene
    drug = ifelse(sapply(strsplit(CARD.Short.Name, "_"), length) == 3,   # Third part (if present): Drug
                  sapply(strsplit(CARD.Short.Name, "_"), `[`, 3), NA)
  )

# View the mutated data
head(aro_index)


# Print the first few rows of each dataframe
print(head(aro_index))
print(head(antibiotics_data))
print(head(pathogens_data))

# Show the column names of each dataframe
print(colnames(aro_index))
print(colnames(antibiotics_data))
print(colnames(pathogens_data))

#Extract pathogen, gene, and drug from 'CARD.Short.Name'
aro_index_clean <- aro_index %>%
  mutate(
    pathogen = sapply(strsplit(CARD.Short.Name, "_"), `[`, 1),  # Extract pathogen
    gene = sapply(strsplit(CARD.Short.Name, "_"), `[`, 2),      # Extract gene
    drug = ifelse(sapply(strsplit(CARD.Short.Name, "_"), length) == 3,   # If present, extract drug
                  sapply(strsplit(CARD.Short.Name, "_"), `[`, 3), NA)
  )

#Merge aro_index_clean with the antibiotics_data and pathogens_data
# For merging with antibiotics_data
merged_data_antibiotics <- left_join(aro_index_clean, antibiotics_data,
                                     by = c("gene" = "AAC.Abbreviation"))

# For merging with pathogens_data
merged_data_pathogens <- left_join(merged_data_antibiotics, pathogens_data,
                                   by = c("pathogen" = "Abbreviation"))

#Remove duplicate rows (if any) and filter out rows where pathogen is empty
cleaned_data <- merged_data_pathogens %>%
  distinct() %>%
  filter(!is.na(pathogen))

#Group by Pathogen, Gene, and Drug, then summarize Antibiotic information
summarized_data <- cleaned_data %>%
  group_by(Pathogen = Pathogen, Gene = gene, Drug = drug) %>%
  summarize(Antibiotic_Info = paste(unique(Molecule), collapse = ", ")) %>%
  arrange(Pathogen, Gene, Drug)

#View the final summarized data
print(head(summarized_data))

# Filter for Klebsiella pneumoniae and CZA(Bug-Drug of Interest)
summarized_data %>%
  filter(Pathogen == "Klebsiella pneumoniae", Drug == "CZA") -> klebsiella_cza_combinations

# View the filtered results
head(klebsiella_cza_combinations)

