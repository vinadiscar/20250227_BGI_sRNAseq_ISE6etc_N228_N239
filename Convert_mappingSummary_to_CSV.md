# Convert mapping summary table (from Bowtie log file) to csv file
- Prepare mapping summary table: SE, RNAseq data (Hlo-infected and control samples) mapped to LGTV or HAZV genome log files will be converted to csv files

## R Script to Extract Mapping Stats and Save to CSV
```sh 
conda activate my_r_env
R

# Set directories
hazv_dir <- "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result/logs"
lgtv_dir <- "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/logs"

# Get all log file paths
hazv_logs <- list.files(hazv_dir, pattern = "_bowtie\\.log$", full.names = TRUE)
lgtv_logs <- list.files(lgtv_dir, pattern = "_bowtie\\.log$", full.names = TRUE)

# Function to extract mapping data from a log file
extract_log_data <- function(file, genome) {
  lines <- readLines(file)
  sample <- gsub("_bowtie\\.log$", "", basename(file))
  
  total_line <- grep("# reads processed:", lines, value = TRUE)
  aligned_line <- grep("# reads with at least one reported alignment:", lines, value = TRUE)

  # Use regular expression to extract the numeric part before any parentheses
  total_reads <- as.numeric(sub(".*: ([0-9]+).*", "\\1", total_line))
  aligned_reads <- as.numeric(sub(".*: ([0-9]+) \\(.*", "\\1", aligned_line))

  percent_aligned <- round((aligned_reads / total_reads) * 100, 6)
  
  data.frame(
    Sample = sample,
    Genome = genome,
    Total_Reads = total_reads,
    Aligned_Reads = aligned_reads,
    Percent_Aligned = percent_aligned,
    stringsAsFactors = FALSE
  )
}

# Extract data
hazv_data <- do.call(rbind, lapply(hazv_logs, extract_log_data, genome = "HAZV"))
lgtv_data <- do.call(rbind, lapply(lgtv_logs, extract_log_data, genome = "LGTV"))

# Combine and save
mapping_summary <- rbind(hazv_data, lgtv_data)
write.csv(mapping_summary, "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/mapping_summary_bowtie1_SE_HlosRNAseq_viralGenome.csv", row.names = FALSE)
```