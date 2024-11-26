# Load required libraries
library(ggplot2)
library(dplyr)

# Set working directory
setwd("C:/Users/ericd/Dropbox/Eric Work/DPA/APHIDS/data/")

# Define the sample ID
sample_id <- "S06"

# Read the depth data from a file
data <- read.table(paste("coverage_",sample_id,".txt",sep=""), header = FALSE, sep = "\t")

# Assign column names to the data
colnames(data) <- c('Contig', 'Position', 'Depth')

# Convert Contig column to a factor to ensure it's treated as a categorical variable
data$Contig <- as.factor(data$Contig)

# Convert Position to kilobases
data$Position <- data$Position / 1000

# Read the list of contigs from a file
contig_list <- read.table('filtered_assemblies_contig_list.txt', header = FALSE, sep = " ")

# Ensure the column names are correct
# If not, rename them to match what your script expects
colnames(contig_list) <- c('SampleID', 'Contig')

# Filter the contig list for the specific sample ID
contigs <- contig_list %>%
  filter(SampleID == sample_id) %>%
  pull(Contig)  # Extract contig names as a vector

# Filter the depth data to include only the specified contigs
filtered_data <- data %>%
  filter(Contig %in% contigs)

# Calculate the median depth for the filtered data
medLine <- median(filtered_data$Depth)

# Read the first BLAST output data
blast_data1 <- read.table('APSE_results.txt', header = FALSE, sep = "\t")

# Assign column names to the first BLAST data
colnames(blast_data1) <- c('SubjectID', 'SubjectStart', 'SubjectEnd', 'QueryID', 'SampleID', 
                           'PercentIdentity', 'AlignmentLength', 'Mismatches', 'GapOpens',
                           'QueryStart', 'QueryEnd', 'SubjectStartPos', 'SubjectEndPos', 'Evalue')

# Filter the first BLAST data for the specific sample ID and contigs
blast_filtered1 <- blast_data1 %>%
  filter(grepl(sample_id, SampleID)) %>%
  mutate(Contig = gsub(".*_>", "", SampleID)) %>%  # Extract contig name from SampleID
  filter(Contig %in% contigs) %>%  # Filter contigs matching the list
  mutate(SubjectStartPos = SubjectStartPos / 1000,  # Convert positions to kilobases
         SubjectEndPos = SubjectEndPos / 1000) %>%
  select(Contig, SubjectStartPos, SubjectEndPos)  # Select relevant columns

# Read the second BLAST output data
blast_data2 <- read.table('PHD5AT_results.txt', header = FALSE, sep = "\t")

# Assign column names to the second BLAST data
colnames(blast_data2) <- c('SubjectID', 'SubjectStart', 'SubjectEnd', 'QueryID', 'SampleID', 
                           'PercentIdentity', 'AlignmentLength', 'Mismatches', 'GapOpens',
                           'QueryStart', 'QueryEnd', 'SubjectStartPos', 'SubjectEndPos', 'Evalue')

# Filter the second BLAST data for the same sample ID and contigs
blast_filtered2 <- blast_data2 %>%
  filter(grepl(sample_id, SampleID)) %>%
  mutate(Contig = gsub(".*_>", "", SampleID)) %>%  # Extract contig name from SampleID
  filter(Contig %in% contigs) %>%  # Filter contigs matching the list
  mutate(SubjectStartPos = SubjectStartPos / 1000,  # Convert positions to kilobases
         SubjectEndPos = SubjectEndPos / 1000) %>%
  select(Contig, SubjectStartPos, SubjectEndPos)  # Select relevant columns

# Read the third BLAST output data
blast_data3 <- read.table('M147_results.txt', header = FALSE, sep = "\t")

# Assign column names to the third BLAST data
colnames(blast_data3) <- c('SubjectID', 'SubjectStart', 'SubjectEnd', 'QueryID', 'SampleID', 
                           'PercentIdentity', 'AlignmentLength', 'Mismatches', 'GapOpens',
                           'QueryStart', 'QueryEnd', 'SubjectStartPos', 'SubjectEndPos', 'Evalue')

# Filter the third BLAST data for the same sample ID and contigs
blast_filtered3 <- blast_data3 %>%
  filter(grepl(sample_id, SampleID)) %>%
  mutate(Contig = gsub(".*_>", "", SampleID)) %>%  # Extract contig name from SampleID
  filter(Contig %in% contigs) %>%  # Filter contigs matching the list
  mutate(SubjectStartPos = SubjectStartPos / 1000,  # Convert positions to kilobases
         SubjectEndPos = SubjectEndPos / 1000) %>%
  select(Contig, SubjectStartPos, SubjectEndPos)  # Select relevant columns

# Read the fourth BLAST output data
blast_data4 <- read.table('P4M47_results.txt', header = FALSE, sep = "\t")

# Assign column names to the fourth BLAST data
colnames(blast_data4) <- c('SubjectID', 'SubjectStart', 'SubjectEnd', 'QueryID', 'SampleID', 
                           'PercentIdentity', 'AlignmentLength', 'Mismatches', 'GapOpens',
                           'QueryStart', 'QueryEnd', 'SubjectStartPos', 'SubjectEndPos', 'Evalue')

# Filter the fourth BLAST data for the same sample ID and contigs
blast_filtered4 <- blast_data4 %>%
  filter(grepl(sample_id, SampleID)) %>%
  mutate(Contig = gsub(".*_>", "", SampleID)) %>%  # Extract contig name from SampleID
  filter(Contig %in% contigs) %>%  # Filter contigs matching the list
  mutate(SubjectStartPos = SubjectStartPos / 1000,  # Convert positions to kilobases
         SubjectEndPos = SubjectEndPos / 1000) %>%
  select(Contig, SubjectStartPos, SubjectEndPos)  # Select relevant columns

# Combine all BLAST filtered data for plotting
blast_filtered1$Type <- 'APSE_results'
blast_filtered2$Type <- 'PHD5AT_results'
blast_filtered3$Type <- 'M147_results'
blast_filtered4$Type <- 'P4M47_results'

# Combine all filtered BLAST data into one data frame for plotting
combined_blast_data <- bind_rows(blast_filtered1, blast_filtered2, blast_filtered3, blast_filtered4)

# Create the ggplot with multiple panels for each specified contig
p <- ggplot(filtered_data, aes(x = Position, y = Depth)) +
  geom_line() +  # Line plot for depth across each position
  geom_hline(yintercept = medLine, col = "red", lty = 2) +  # Add a horizontal line for the median depth
  geom_segment(data = combined_blast_data, aes(x = SubjectStartPos, xend = SubjectEndPos, y = 0, yend = 0, color = Type), size = 3) +  # Add BLAST hits with different colors
  facet_wrap(~ Contig, scales = 'free_x') +  # Create separate panels for each contig
  theme_minimal() +  # Apply a minimal theme
  labs(title = "Read Mapping Depth and BLAST Hits Across Selected Contigs",
       x = "Position (Kb)",  # Update x-axis label to reflect Kb
       y = "Depth (log10 scale)") +  # Add labels
  scale_y_log10() +  # Change y-axis to log10 scale
  scale_color_manual(values = c("blue", "green", "orange", "red"), name = "BLAST Hit Type")  # Add a legend for BLAST hit types

# Save the plot to a single file
ggsave(filename = paste("diagnostics_",sample_id,".png",sep=""), plot = p, width = 25, height = 15, units = "cm", bg = "white")
