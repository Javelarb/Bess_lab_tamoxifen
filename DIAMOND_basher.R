# Load necessary libraries
library(tidyverse)
library(readr)
library(ggpubr)

setwd("/Users/juli.avelar-barragan/Desktop/Project_bioinformatics/Bess_lab_tamoxifen/")

# Get a list of all txt files in the directory
file_list <- list.files(path = "HMGI3013_diamond_out/", pattern = "\\.txt$", full.names = TRUE)

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

#write.table(Gus_RC_wide, "HMGC279_RPM_counts.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Supplemental coverage figure
cov_df = combined_df %>%
  filter(X2 == "SRS024549.1421-T1-C", 
         !file_name %in% c("2126", "22b", "2621", "Stan"))

gene_length = max(cov_df$X10)

calculate_sample_coverage <- function(sample_df, gene_length) {
  sample_df %>%
    rowwise() %>%
    do(data.frame(position = seq(.$X9, .$X10))) %>%
    ungroup() %>%
    count(position) %>%
    complete(position = 1:gene_length, fill = list(n = 0)) %>%
    select(position, coverage = n)
}

calculate_all_samples_coverage <- function(df) {
  df %>%
    group_by(file_name) %>%
    group_modify(~ calculate_sample_coverage(.x, gene_length)) %>%
    ungroup()
}

all_samples_coverage <- calculate_all_samples_coverage(cov_df)

# Calculate average coverage across all samples
average_coverage <- all_samples_coverage %>%
  group_by(position) %>%
  summarize(avg_coverage = mean(coverage))

# Calculate overall average coverage
overall_avg_coverage <- mean(average_coverage$avg_coverage)

# Create the plot
A = ggplot(average_coverage, aes(x = position, y = avg_coverage)) +
  geom_line() +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("Overall Gene Coverage: %.2f", overall_avg_coverage),
           hjust = 1, vjust = 1) +
  labs(x = "Position in Gene", y = "Average Coverage", 
       title = "B. fragilis", subtitle = "SRS024549.1421-T1-C") +
  theme_bw()

cov_df2 = combined_df %>%
  filter(X2 == "SRS024388.120461-T1-C", 
         !file_name %in% c("2126", "22b", "2621", "Stan"))

all_samples_coverage2 <- calculate_all_samples_coverage(cov_df2)

# Calculate average coverage across all samples
average_coverage2 <- all_samples_coverage2 %>%
  group_by(position) %>%
  summarize(avg_coverage = mean(coverage))

# Calculate overall average coverage
overall_avg_coverage2 <- mean(average_coverage2$avg_coverage)

# Create the plot
B = ggplot(average_coverage2, aes(x = position, y = avg_coverage)) +
  geom_line() +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("Overall Gene Coverage: %.2f", overall_avg_coverage2),
           hjust = 1, vjust = 1) +
  labs(x = "Position in Gene", y = "Average Coverage", 
       title = "B. fragilis", subtitle = "SRS024388.120461-T1-C") +
  theme_bw()

ggsave(filename = "Supplemental_coverage_plot.png", plot = ggarrange(A, B, ncol = 1, labels = c("A", "B")),
      dpi = 600, bg = "white", width = 6.54, height = 5)
