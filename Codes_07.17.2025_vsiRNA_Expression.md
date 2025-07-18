# Analyze vsiRNA expression
### Goals:
###### 1. Size distribution of vsiRNAs 
###### 2. Positional distribution on the viral genome 

## Part 1. Map Reads to viral genome (HAZV & LGTV)
- workflow (map to HAZV genome) can be found in [github@vinadiscar "Map small RNAseq data Hlo to HAZV genome using Bowtie 1"](https://github.com/vinadiscar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/blob/main/BowtieMap_SE_sRNAseqHlo_to_HAZVgenome_update2.md)
- workflow (map to LGTV genome) can be found in [github@vinadiscar "Map small RNAseq data Hlo to HAZV genome using Bowtie 1"](https://github.com/vinadiscar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/blob/main/Bowtie_map_sRMAseqHlo_to_LGTVgenome_update.md)
- Output files: SAM files

- next step: Convert Sam file to BAM files

#### Convert SAM files to BAM files

- ###### WRITE SCRIPT (convert_sam_to_bam_HAZV.sh)
```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update
vi convert_sam_to_bam_HAZV.sh
```
- convert_sam_to_bam_HAZV.sh
```sh
#!/bin/bash

# Set working directory
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update

# Loop through samples N228 to N239
for i in {228..239}; do
    echo "Processing sample N$i..."

    # Convert SAM to BAM
    samtools view -bS HAZVgenome_N${i}_L1.sam > HAZVgenome_N${i}_L1.bam

    # Sort BAM
    samtools sort -o HAZVgenome_N${i}_L1.sorted.bam HAZVgenome_N${i}_L1.bam

    # Index sorted BAM
    samtools index HAZVgenome_N${i}_L1.sorted.bam

    echo "Sample N$i done."
done

echo "All SAM files processed."
```
- ###### RUN: convert_sam_to_bam_HAZV.sh

```sh
chmod +x convert_sam_to_bam_HAZV.sh
./convert_sam_to_bam_HAZV.sh
```
- ##### ***Do the same with LGTV Sam files (convert_sam_to_bam_LGTV.sh)


- ###### WRITE SCRIPT (convert_sam_to_bam_LGTV.sh)
```sh

cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update
vi convert_sam_to_bam_LGTV.sh
```
- convert_sam_to_bam_LGTV.sh
```sh
#!/bin/bash

# Set working directory
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update

# Loop through samples N228 to N239
for i in {228..239}; do
    echo "Processing sample N$i..."

    # Convert SAM to BAM
    samtools view -bS LGTVgenome_N${i}_L1.sam > LGTVgenome_N${i}_L1.bam

    # Sort BAM
    samtools sort -o LGTVgenome_N${i}_L1.sorted.bam LGTVgenome_N${i}_L1.bam

    # Index sorted BAM
    samtools index LGTVgenome_N${i}_L1.sorted.bam

    echo "Sample N$i done."
done

echo "All SAM files processed."
```
- ###### RUN: convert_sam_to_bam_LGTV.sh

```sh
chmod +x convert_sam_to_bam_LGTV.sh
./convert_sam_to_bam_LGTV.sh
```

## Part 2 Size Distribution of vsiRNAs and Genomic Distribution
###### Size Distribution
- read lengths per segment
###### Genomic Distribution
- Read distribution along each viral genome segment
- mapping the read coverage along each viral segment


```sh
R
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)
library(dplyr)
library(tidyr)
```
### A. HAZV 
```sh

setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update")

# 1. Define BAM File Paths
bam_dir_HAZV <- "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update" # set directory

# List BAM files
bam_files_HAZV <- list.files(bam_dir_HAZV, pattern = "\\.sorted\\.bam$", full.names = TRUE)

# 2. Extract Read Length and Position Information
extract_vsiRNA_info_HAZV <- function(bam_file_HAZV) {
  aln <- readGAlignments(bam_file_HAZV, use.names = TRUE)
  
  df_HAZV <- data.frame(
    seqnames = as.character(seqnames(aln)),       # segment name
    start = start(aln),
    width = width(aln)                            # read length
  )
  df_HAZV$sample <- basename(bam_file_HAZV)
  return(df_HAZV)
}

# Apply to all BAM files
vsiRNA_df_HAZV <- do.call(rbind, lapply(bam_files_HAZV, extract_vsiRNA_info_HAZV))

write.csv(vsiRNA_df_HAZV, "vsiRNA_all_reads_HAZV.csv", row.names = FALSE) # Basic Table of All Mapped Reads


# 3. Plot Histogram of Read Sizes per Segment
# Plot for all segments and samples
pdf("vsiRNA_read_length_histogram_HAZV.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_HAZV, aes(x = width)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ seqnames, scales = "free_y") +
  labs(title = "vsiRNA Read Length Distribution per Segment_HAZV",
       x = "Read Length (nt)", y = "Count") +
  theme_minimal(base_size = 14)
dev.off()

# or when log-scaling is applied

pdf("vsiRNA_read_length_histogram_HAZV_logY.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_HAZV, aes(x = width)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  scale_y_log10() +  # <- Log scale added here
  facet_wrap(~ seqnames, scales = "free_y") +
  labs(title = "vsiRNA Read Length Distribution per Segment_HAZV (log-scaled)",
       x = "Read Length (nt)", y = "Log10(Read Count)") +
  theme_minimal(base_size = 14)
dev.off()


#  4. Plot Distribution Along Viral Genome
# Distribution of vsiRNAs across segments

pdf("vsiRNA_distribution_on_genome_HAZV.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_HAZV, aes(x = start)) +
  geom_histogram(binwidth = 100, fill = "tomato", color = "black") +
  facet_wrap(~ seqnames, scales = "free_x") +
  labs(title = "vsiRNA Mapping Distribution on HAZV Viral Genome",
       x = "Genomic Position", y = "Read Count") +
  theme_minimal(base_size = 14)
dev.off()

# Save Summary Key information of vsiRNA in CSV file

vsiRNA_summary_HAZV <- vsiRNA_df_HAZV %>%
  group_by(sample, seqnames) %>%
  summarise(
    total_reads = n(),
    mean_read_length = mean(width),
    most_common_length = as.numeric(names(sort(table(width), decreasing = TRUE)[1])),
    min_length = min(width),
    max_length = max(width)
  ) %>%
  arrange(sample, seqnames)

write.csv(vsiRNA_summary_HAZV, "vsiRNA_summary_per_sample_and_segment_HAZV.csv", row.names = FALSE)

```

### B. LGTV

```sh

# 1. Define BAM File Paths
bam_dir_LGTV <- "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update" # set directory

# List BAM files
bam_files_LGTV <- list.files(bam_dir_LGTV, pattern = "\\.sorted\\.bam$", full.names = TRUE)

# 2. Extract Read Length and Position Information
extract_vsiRNA_info_LGTV <- function(bam_file_LGTV) {
  aln <- readGAlignments(bam_file_LGTV, use.names = TRUE)
  
  df_LGTV <- data.frame(
    seqnames = as.character(seqnames(aln)),       # segment name
    start = start(aln),
    width = width(aln)                            # read length
  )
  df_LGTV$sample <- basename(bam_file_LGTV)
  return(df_LGTV)
}

# Apply to all BAM files
vsiRNA_df_LGTV <- do.call(rbind, lapply(bam_files_LGTV, extract_vsiRNA_info_LGTV))

write.csv(vsiRNA_df_LGTV, "vsiRNA_all_reads_LGTV.csv", row.names = FALSE) # Basic Table of All Mapped Reads


# 3. Plot Histogram of Read Sizes per Segment
# Plot for all segments and samples
pdf("vsiRNA_read_length_histogram_LGTV.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_LGTV, aes(x = width)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ seqnames, scales = "free_y") +
  labs(title = "vsiRNA Read Length Distribution per Segment_LGTV",
       x = "Read Length (nt)", y = "Count") +
  theme_minimal(base_size = 14)
dev.off()

# or when log-scaling is applied

pdf("vsiRNA_read_length_histogram_LGTV_logY.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_LGTV, aes(x = width)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  scale_y_log10() +  # <- Log scale added here
  facet_wrap(~ seqnames, scales = "free_y") +
  labs(title = "vsiRNA Read Length Distribution per Segment_LGTV",
       x = "Read Length (nt)", y = "Log10(Read Count)") +
  theme_minimal(base_size = 14)
dev.off()

#  4. Plot Distribution Along Viral Genome
# Distribution of vsiRNAs across segments

pdf("vsiRNA_distribution_on_genome_LGTV.pdf", width = 10, height = 6)
ggplot(vsiRNA_df_LGTV, aes(x = start)) +
  geom_histogram(binwidth = 100, fill = "tomato", color = "black") +
  facet_wrap(~ seqnames, scales = "free_x") +
  labs(title = "vsiRNA Mapping Distribution on LGTV Viral Genome",
       x = "Genomic Position", y = "Read Count") +
  theme_minimal(base_size = 14)
dev.off()

# Save Summary Key information of vsiRNA in CSV file

vsiRNA_summary_LGTV <- vsiRNA_df_LGTV %>%
  group_by(sample, seqnames) %>%
  summarise(
    total_reads = n(),
    mean_read_length = mean(width),
    most_common_length = as.numeric(names(sort(table(width), decreasing = TRUE)[1])),
    min_length = min(width),
    max_length = max(width)
  ) %>%
  arrange(sample, seqnames)

write.csv(vsiRNA_summary_LGTV, "vsiRNA_summary_LGTV.csv", row.names = FALSE)
```







