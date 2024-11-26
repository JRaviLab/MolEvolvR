# Load required libraries
library(dplyr)
library(readr)

# File paths
original_file <- "CARD_data/aro_index.tsv"
cleaned_file <- "CARD_data/resistance_profile_data.tsv"

# Load datasets
aro_index <- read_delim("CARD_data/aro_index.tsv", delim = "\t", col_names = TRUE)
resistance_profile_data <- read_delim("CARD_data/resistance_profile_data.tsv", delim = "\t", col_names = TRUE)

# Extract the first 10rows from each file
aro_index_snippet <- head(aro_index, 10)
resistance_profile_data_snippet <- head(resistance_profile_data, 10)

# View the pre-cleanup snippet
View(aro_index_snippet)

# View the post-cleanup snippet
View(resistance_profile_data_snippet)

