setwd("C:\\Users\\Administrator\\Desktop\\card-data SNPs")
snps_data <- read.delim("snps.txt", header = TRUE, sep = "\t")
head(snps_data)
# Load tibble package
library(tibble)
# Read the SNP data into a data frame
snps_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert data frame to tibble
snps_tibble <- as_tibble(snps_data)

# View the tibble
print(snps_tibble)

# Load the dplyr package
library(dplyr)

# Define the pattern to search for in the Name column
search_pattern <- "Klebsiella pneumoniae lamB with mutations conferring resistance to ceftazidime-avibactam"

# Filter the tibble for entries that match the search pattern in the Name column
filtered_snps <- snps_tibble %>%
  filter(grepl(search_pattern, Name))

# View the filtered result
print(filtered_snps)

# Save the filtered tibble to a CSV file
write.csv(filtered_snps, file = "Klebsiella_pneumoniae_lamB_resistance.csv", row.names = FALSE)


