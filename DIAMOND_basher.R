# Load necessary libraries
library(tidyverse)
library(readr)

setwd("/Users/juli.avelar-barragan/Desktop/Project_bioinformatics/Bess_lab_tamoxifen/")

# Get a list of all txt files in the directory
file_list <- list.files(path = "HMGC279_diamond_out/", pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
df_list <- list()

# Loop through each file and read it into a data frame
for (file in file_list) {
  # Read the file into a data frame
  df <- read_delim(file, delim = "\t", col_names = F, comment = "#")
  
  # Add a column with the file name
  df <- df %>% mutate(file_name = sub("\\.txt$", "", basename(file)))
  
  # Add the data frame to the list
  df_list[[length(df_list) + 1]] <- df
}

# Combine all data frames into a single data frame
combined_df <- bind_rows(df_list)

#Count the number of reads. 
#Reads which map to multiple features will only have the top feature counted (highest %ID and lowest e-value)
#If there is still multiple features, then the read count is divided equally.

rc_info = read.csv("Readinfo.csv") %>%
  mutate(RPM = Host_removed/1e6) %>%
  select(SampleID, RPM)

GUS_rc <- combined_df %>%
  group_by(file_name, X1) %>%
  filter(X3 == max(X3) & X11 == min(X11) & X12 == max(X12)) %>%
  mutate(Count = 1 / n()) %>%
  ungroup() %>%
  group_by(file_name, X2) %>%
  summarize(FeatureCount = sum(Count), .groups = 'drop') %>%
  ungroup() %>%
  group_by(file_name) %>%
  mutate(TotalCount = sum(FeatureCount)) %>%
  right_join(rc_info, by = join_by("file_name" == "SampleID")) %>%
  mutate(NormFeatureCount = FeatureCount/RPM)

Gus_RC_wide = GUS_rc %>%
  pivot_wider(id_cols = file_name, names_from = X2, values_from = NormFeatureCount) %>%
  column_to_rownames(var = "file_name") %>%
  mutate_all(~ coalesce(., 0)) %>%
  rownames_to_column(var = "SampleID")

write.table(Gus_RC_wide, "HMGC279_RPM_counts.txt", quote = F, sep = "\t", row.names = F, col.names = T)
