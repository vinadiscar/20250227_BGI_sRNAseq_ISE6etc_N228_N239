# Strand-specific vsiRNA size distribution: 5'nucleotide Analysis

FASTQ
  ↓
Mapping to HAZV genome (Bowtie1)
  ↓
SAM file
  ↓
Count vsiRNAs (awk/samtools)
  ↓
vsi_counts_all.tsv   ← raw counts
  ↓
CPM normalization (CPM is calculated based on the number of reads mapped to the host genome)
  ↓
Plots/figures

#### Feb 5, 2026 

## II HAZV-Infected (Updated): Each Day/Time Point Update Version 2 
- update CPM normalization method (CPM is calculated based on the number of reads mapped to the host genome)



### 1. Generate vsiRNA table from viral mapping
```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints

# Map sample IDs to Day
declare -A sam_day
sam_day=(
  [N230]="Day1"
  [N233]="Day2"
  [N236]="Day3"
  [N239]="Day4"
)

# Create single raw file
> HAZV_vsi_raw_all_days.tsv  # empty file to start

for sam in "${!sam_day[@]}"; do
    day=${sam_day[$sam]}
    echo "Processing $day ($sam)..."

    samtools view -F 4 "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_${sam}_L1.sam" | \
    awk -v day="$day" '{
        ref=$3; seq=$10; len=length(seq);      
        strand = (and($2,16)) ? "-" : "+";
        print day "\t" ref "\t" seq "\t" len "\t" strand
    }' >> HAZV_vsi_raw_all_days.tsv
done

# Collapse identical sequences per day and segment
LC_ALL=C sort HAZV_vsi_raw_all_days.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1}' \
> HAZV_vsi_counts_all_days.tsv

```
- Output files:
    - HAZV_vsi_raw_all_days.tsv
    - HAZV_vsi_counts_all_days.tsv


### 2. Get host-mapped read counts
```sh
conda activate samtoools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints_update

echo -e "Day\thost_mapped_reads" > HAZV_host_mapped_reads_update.tsv

samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N230_L1Aligned.sortedByCoord.out.bam | awk '{print "Day1\t"$1}' >> HAZV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N233_L1Aligned.sortedByCoord.out.bam | awk '{print "Day2\t"$1}' >> HAZV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N236_L1Aligned.sortedByCoord.out.bam | awk '{print "Day3\t"$1}' >> HAZV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N239_L1Aligned.sortedByCoord.out.bam | awk '{print "Day4\t"$1}' >> HAZV_host_mapped_reads_update.tsv
```

### 3. R analysis: 5' nucleotide and host-normalized CPM
```sh
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints_update")

# Read vsiRNA counts (viral mapping)
vsi <- read_tsv(
  "HAZV_vsi_counts_all_days.tsv",
  col_names = c("Day","segment","seq","length","strand","count")
)

# Read host-mapped read counts
host_counts <- read_tsv("HAZV_host_mapped_reads_update.tsv")

head(host_counts)
# Day    host_mapped_reads
#  <chr>             <dbl>
# 1 Day1          137431056
# 2 Day2          174867073
# 3 Day3          138264769
# 4 Day4          155666927

# Filter length & extract correct 5' nucleotide
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = case_when(
      strand == "+" ~ substr(seq, 1, 1),
      strand == "-" ~ substr(seq, nchar(seq), nchar(seq))
    ),
    five_prime = toupper(five_prime),
    five_prime = recode(five_prime, T = "U"),
    five_prime = factor(five_prime, levels = c("U","A","C","G"))
  ) %>%
  filter(!is.na(five_prime))

# ---- CPM calculations ----

# Host-normalized CPM
vsi <- vsi %>%
  left_join(host_counts, by = "Day") %>%
  mutate(
    CPM_host = (count / host_mapped_reads) * 1e6
  )
head(vsi)

write_tsv(vsi, "HAZV_vsiRNA_CPM_host_update.tsv")

# Summarize (use CPM_host)
vsi_summary <- vsi %>%
  group_by(Day, segment, length, strand, five_prime) %>%
  summarise(
    CPM_host = sum(CPM_host),
    .groups = "drop"
  )
head (vsi_summary)
write_tsv(vsi_summary, "HAZV_vsiRNA_CPM_host_summary_update.tsv")

# Mirror plot values
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand == "-", -CPM_host, CPM_host))
head(vsi_mirror)

write_tsv(vsi_mirror, "HAZV_vsiRNA_CPM_mirror_host_update.tsv")
```

### 4. Generate Plots

```sh
# Use host-normalized mirror data
plot_data <- vsi_mirror

days <- unique(plot_data$Day)

for (day in days) {

  segments <- unique(filter(plot_data, Day == day)$segment)

  for (seg in segments) {

    p <- ggplot(
      filter(plot_data, Day == day, segment == seg),
      aes(x = length, y = CPM_mirror, fill = five_prime)
    ) +
      geom_col(width = 0.8) +
      scale_x_continuous(
        breaks = 18:30,
        limits = c(17.5, 30.5)
      ) +
      scale_y_continuous(
        labels = function(x) comma(abs(x)),
        name = "CPM (normalized to host genome)"
      ) +
      scale_fill_manual(
        values = c(
          "U" = "#E41A1C",
          "A" = "#4DAF4A",
          "C" = "#377EB8",
          "G" = "#984EA3"
        )
      ) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      labs(
        title = paste("HAZV vsiRNAs –", seg, day),
        x = "vsiRNA length (nt)",
        fill = "5′ nucleotide"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )

    ggsave(
      filename = paste0("HAZV_vsiRNA_", seg, "_hostCPM_", day, ".pdf"),
      plot = p,
      width = 6,
      height = 8,
      device = cairo_pdf
    )
  }
}

```

February 5, 2026

## II LGTV-Infected (Update): Each Day/Time Point Update Version 2 

- update CPM normalization method (CPM is calculated based on the number of reads mapped to the host genome)

### 1. Generate vsiRNA table from viral mapping

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints_update

> LGTV_vsi_raw_update.tsv

declare -A sam_day=(
  [N230]="Day1"
  [N233]="Day2"
  [N236]="Day3"
  [N239]="Day4"
)

for sam in "${!sam_day[@]}"; do
    day=${sam_day[$sam]}
    echo "Processing $day ($sam)..."

    samtools view -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_${sam}_L1.sam | \
    awk -v day="$day" '{
        seq=$10;
        len=length(seq);
        strand = (and($2,16)) ? "-" : "+";
        print day "\t" seq "\t" len "\t" strand
    }' >> LGTV_vsi_raw_update.tsv
done

# Collapse identical reads per day
LC_ALL=C sort LGTV_vsi_raw_update.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' \
> LGTV_vsi_counts_update.tsv

```

### 2. Get host mapped read counts

```sh 
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints_update

echo -e "Day\thost_mapped_reads" > LGTV_host_mapped_reads_update.tsv

samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N229_L1Aligned.sortedByCoord.out.bam | awk '{print "Day1\t"$1}' >> LGTV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N232_L1Aligned.sortedByCoord.out.bam | awk '{print "Day2\t"$1}' >> LGTV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N235_L1Aligned.sortedByCoord.out.bam | awk '{print "Day3\t"$1}' >> LGTV_host_mapped_reads_update.tsv
samtools view -c -F 4 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR_update/HaeL_N238_L1Aligned.sortedByCoord.out.bam | awk '{print "Day4\t"$1}' >> LGTV_host_mapped_reads_update.tsv

```
### 3. R analysis: 5' nucleotide and host-normalized CPM

```sh
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints_update")

# LGTV vsiRNA counts (viral mapping)
vsi <- read_tsv(
  "LGTV_vsi_counts_update.tsv",
  col_names = c("Day","seq","length","strand","count")
)

# Host-mapped read counts
host_counts <- read_tsv("LGTV_host_mapped_reads_update.tsv")

head(host_counts)
# A tibble: 4 × 2
#  Day   host_mapped_reads
# <chr>             <dbl>
# Day1          140368650
# Day2          208801730
# Day3          114972138
# Day4          128381876

# Filter lengths & extract correct 5' nucleotide
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = case_when(
      strand == "+" ~ substr(seq, 1, 1),
      strand == "-" ~ substr(seq, nchar(seq), nchar(seq))
    ),
    five_prime = toupper(five_prime),
    five_prime = recode(five_prime, T = "U"),
    five_prime = factor(five_prime, levels = c("U","A","C","G"))
  ) %>%
  filter(!is.na(five_prime))
head(vsi)

# Host-normalized CPM
vsi <- vsi %>%
  left_join(host_counts, by = "Day") %>%
  mutate(
    CPM_host = (count / host_mapped_reads) * 1e6
  )
head(vsi)

write_tsv(vsi, "LGTV_vsiRNA_CPM_host_update.tsv")

# Summarize for Plotting

vsi_summary <- vsi %>%
  group_by(Day, length, strand, five_prime) %>%
  summarise(
    CPM_host = sum(CPM_host),
    .groups = "drop"
  )

write_tsv(vsi_summary, "LGTV_vsiRNA_CPM_summary_host_update.tsv")

# Mirror values
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand == "-", -CPM_host, CPM_host))

write_tsv(vsi_mirror, "LGTV_vsiRNA_CPM_mirror_host_update.tsv")
```
### 4. Generate plot

```sh
days <- unique(vsi_mirror$Day)

for (day in days) {

  p <- ggplot(
    filter(vsi_mirror, Day == day),
    aes(x = length, y = CPM_mirror, fill = five_prime)
  ) +
    geom_col(width = 0.8) +
    scale_x_continuous(
      breaks = 18:30,
      limits = c(17.5, 30.5)
    ) +
    scale_y_continuous(
      labels = function(x) comma(abs(x)),
      name = "CPM (normalized to host genome)"
    ) +
    scale_fill_manual(
      values = c(
        "U" = "#E41A1C",
        "A" = "#4DAF4A",
        "C" = "#377EB8",
        "G" = "#984EA3"
      )
    ) +
    geom_hline(yintercept = 0) +
    labs(
      title = paste("LGTV vsiRNAs", day),
      x = "vsiRNA length (nt)",
      fill = "5′ nucleotide"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(
    filename = paste0("lgtv_vsiRNA_hostCPM_", day, ".pdf"),
    plot = p,
    width = 6,
    height = 8,
    device = cairo_pdf
  )
}
```

### Previous Code (not valid)
## II HAZV-Infected (Update): Each Day/Time Point
#### 1. Generate Input Frame (Extract Raw Reads)

```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints

# Map sample IDs to Day
declare -A sam_day
sam_day=(
  [N230]="Day1"
  [N233]="Day2"
  [N236]="Day3"
  [N239]="Day4"
)

# Create single raw file
> HAZV_vsi_raw_all_days.tsv  # empty file to start

for sam in "${!sam_day[@]}"; do
    day=${sam_day[$sam]}
    echo "Processing $day ($sam)..."

    samtools view -F 4 "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_${sam}_L1.sam" | \
    awk -v day="$day" '{
        ref=$3; seq=$10; len=length(seq);      
        strand = (and($2,16)) ? "-" : "+";
        print day "\t" ref "\t" seq "\t" len "\t" strand
    }' >> HAZV_vsi_raw_all_days.tsv
done

# Collapse identical sequences per day and segment
LC_ALL=C sort HAZV_vsi_raw_all_days.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1}' \
> HAZV_vsi_counts_all_days.tsv

```

### 2. Generate Plots Per Day

```sh
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints")

# Read combined counts
vsi <- read_tsv("HAZV_vsi_counts_all_days.tsv",
                col_names = c("Day","segment","seq","length","strand","count"))
head(vsi)

# Filter lengths & extract 5' nucleotide
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = toupper(substr(seq,1,1)),
    five_prime = recode(five_prime, T="U"),
    five_prime = factor(five_prime, levels=c("U","A","C","G"))
  ) %>%
  filter(!is.na(five_prime))
head(vsi)

# CPM per Day & Segment
vsi <- vsi %>%
  group_by(Day, segment) %>%
  mutate(
    total_mapped_reads = sum(count),
    CPM = (count / total_mapped_reads) * 1e6
  ) %>%
  ungroup()
head(vsi)

write_tsv(vsi, "HAZV_vsiRNA_CPM_per_segment_all_days.tsv")

# Summarize CPM per Day, Segment, Length, Strand, 5' nucleotide
vsi_summary <- vsi %>%
  group_by(Day, segment, length, strand, five_prime) %>%
  summarise(CPM = sum(CPM), .groups="drop")
head(vsi_summary)
write_tsv(vsi_summary, "HAZV_vsiRNA_CPM_summary_all_days.tsv")

# Mirror values for minus strand
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand=="-", -CPM, CPM))

write_tsv(vsi_mirror, "HAZV_vsiRNA_CPM_mirror_all_days.tsv")

# Generate mirror plots per Day & Segment
days <- unique(vsi_mirror$Day)
for (day in days) {
  segments <- unique(filter(vsi_mirror, Day==day)$segment)
  for (seg in segments) {
    p <- ggplot(filter(vsi_mirror, Day==day, segment==seg),
                aes(x=length, y=CPM_mirror, fill=five_prime)) +
      geom_col(width=0.8) +
      scale_x_continuous(breaks=18:30) +
      scale_y_continuous(labels=function(x) scales::comma(abs(x)), name="CPM") +
      scale_fill_manual(values=c("U"="#E41A1C","A"="#4DAF4A",
                                 "C"="#377EB8","G"="#984EA3")) +
      labs(title=paste("HAZV segment", seg, day),
           x="vsiRNA length (nt)", fill="5' nucleotide") +
      theme_classic() +
      geom_hline(yintercept=0)

    ggsave(filename=paste0("HAZV_vsiRNA_", seg, "_mirror_", day, ".pdf"),
           plot=p, width=6, height=8, device=cairo_pdf)
  }
}

```
## IA. HAZV-infected samples: All Days (Updated)
- include HAZV-infected samples only

#### 1. Generate Input Dataframe (Extract raw reads)

```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific

# HAZV-infected samples only
for sam in \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N230_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N233_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N236_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N239_L1.sam
do
  samtools view -F 4 "$sam" | \
  awk '{
    ref=$3;               # segment ID
    seq=$10;              # read sequence
    len=length(seq);      
    strand = (and($2,16)) ? "-" : "+";
    print ref "\t" seq "\t" len "\t" strand
  }'
done > HAZV_vsi_raw_all_update2.tsv

# Collapse identical reads (unique sequences per segment)
LC_ALL=C sort HAZV_vsi_raw_all_update2.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' \
> HAZV_vsi_counts_all_update2.tsv

```
##### Output files:
- HAZV_vsi_raw_all_update2.tsv = This is raw per-read information extracted from the SAM files. 
- HAZV_vsi_counts_all_update2.tsv # this .tsv file is collapsed counts after sort | uniq -c. Each row = unique small RNA species, with how many times it appeared. It will be used to create weighted histogram

#### 2. Generate specific-stranded vsiRNA distribution Histogram

```sh
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

# Set working directory
setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific")

# Read collapsed counts
vsi <- read_tsv(
  "HAZV_vsi_counts_all_update2.tsv",
  col_names = c("segment", "seq", "length", "strand", "count")
)
head(vsi)

# Filter for canonical vsiRNA lengths & extract 5' nucleotide
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = toupper(substr(seq, 1, 1)),
    five_prime = recode(five_prime, T = "U"),
    five_prime = factor(five_prime, levels = c("U", "A", "C", "G"))
  ) %>%
  filter(!is.na(five_prime))
head(vsi)

# CPM normalization per segment
vsi <- vsi %>%
  group_by(segment) %>%
  mutate(
    total_mapped_reads = sum(count),
    CPM = (count / total_mapped_reads) * 1e6
  ) %>%
  ungroup()
head(vsi)

write_tsv(vsi, "HAZV_vsiRNA_CPM_per_segment_update2.tsv")

# Summarize CPM per segment, length, strand, 5′ nucleotide
vsi_summary <- vsi %>%
  group_by(segment, length, strand, five_prime) %>%
  summarise(CPM = sum(CPM), .groups = "drop")
head(vsi_summary)

write_tsv(vsi_summary, "HAZV_vsiRNA_CPM_summary_update2.tsv")

# Mirror plot values for minus strand
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand == "-", -CPM, CPM))
head(vsi_mirror)

write_tsv(vsi_mirror, "HAZV_vsiRNA_CPM_mirror_update2.tsv")

# Generate mirror plots per segment
segments <- unique(vsi_mirror$segment)

for (seg in segments) {

  p <- ggplot(
    filter(vsi_mirror, segment == seg),
    aes(x = length, y = CPM_mirror, fill = five_prime)
  ) +
    geom_col(width = 0.8) +
    scale_x_continuous(breaks = 18:30) +
    scale_y_continuous(
      labels = function(x) scales::comma(abs(x)),
      name = "CPM"
    ) +
    scale_fill_manual(
      values = c("U" = "#E41A1C", "A" = "#4DAF4A",
                 "C" = "#377EB8", "G" = "#984EA3")
    ) +
    labs(
      title = paste("HAZV segment", seg),
      x = "vsiRNA length (nt)",
      fill = "5' nucleotide"
    ) +
    theme_classic() +
    geom_hline(yintercept = 0)

  ggsave(
    filename = paste0("HAZV_vsiRNA_", seg, "_mirror_update2.pdf"),
    plot = p,
    width = 6,
    height = 8,
    device = cairo_pdf
  )
}
```

# II LGTV-Infected (Update): Each Day/Time Point

```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints

# Map sample IDs to Days
declare -A sam_day
sam_day=(
  [N229]="Day1"
  [N232]="Day2"
  [N235]="Day3"
  [N238]="Day4"
)

# Create one combined raw file
> LGTV_vsi_raw_all_days.tsv  # empty file

for sam in "${!sam_day[@]}"; do
    day=${sam_day[$sam]}
    echo "Processing $day ($sam)..."

    samtools view -F 4 "/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_${sam}_L1.sam" | \
    awk -v day="$day" '{
        seq=$10;
        len=length(seq);
        strand = (and($2,16)) ? "-" : "+";
        print day "\t" seq "\t" len "\t" strand
    }' >> LGTV_vsi_raw_all_days.tsv
done

# Collapse counts (unique sequences) for all days
LC_ALL=C sort LGTV_vsi_raw_all_days.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' > LGTV_vsi_counts_all_days.tsv

```

## 2. Generate Plots

```sh
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints")

# Read counts
vsi <- read_tsv("LGTV_vsi_counts_all_days_fixed.tsv",
                col_names=c("Day","seq","length","strand","count"))
head(vsi)

# Filter lengths & extract 5' nucleotide
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = toupper(substr(seq,1,1)),
    five_prime = recode(five_prime, T="U"),
    five_prime = factor(five_prime, levels=c("U","A","C","G"))
  ) %>%
  filter(!is.na(five_prime))

# CPM per Day
vsi <- vsi %>%
  group_by(Day) %>%
  mutate(
    total_mapped_reads = sum(count),
    CPM = (count / total_mapped_reads) * 1e6
  ) %>%
  ungroup()
head(vsi)

write_tsv(vsi, "LGTV_vsiRNA_CPM_per_sequence_all_days.tsv")

# Summarize CPM per Day / Length / Strand / 5′ nucleotide
vsi_summary <- vsi %>%
  group_by(Day, length, strand, five_prime) %>%
  summarise(CPM = sum(CPM), .groups="drop")

write_tsv(vsi_summary, "LGTV_vsiRNA_CPM_summary_all_days.tsv")

# Mirror values for minus strand
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand=="-", -CPM, CPM))

write_tsv(vsi_mirror, "LGTV_vsiRNA_CPM_mirror_all_days.tsv")

# Generate mirror plots per Day
days <- unique(vsi_mirror$Day)
for (day in days) {
  p <- ggplot(filter(vsi_mirror, Day==day),
              aes(x=length, y=CPM_mirror, fill=five_prime)) +
    geom_col(width=0.8) +
    scale_x_continuous(breaks=18:30) +
    scale_y_continuous(labels=function(x) comma(abs(x)), name="CPM") +
    scale_fill_manual(values=c("U"="#E41A1C","A"="#4DAF4A",
                               "C"="#377EB8","G"="#984EA3")) +
    labs(title=paste("LGTV vsiRNA distribution", day),
         x="vsiRNA length (nt)", fill="5' nucleotide") +
    theme_classic() +
    geom_hline(yintercept=0)

  ggsave(filename=paste0("LGTV_vsiRNA_mirror_", day, ".pdf"),
         plot=p, width=6, height=8, device=cairo_pdf)
}

```
## II. LGTV-infected samples (Update): All Days
- use only LGTV-infected samples

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific
conda activate samtools_env

# Extract sequence, strand, and length in .sam files
for sam in \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_N229_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_N232_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_N235_L1.sam \
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_LGTVgenome_bowtie1/mapping_result_update/LGTVgenome_N238_L1.sam
do
  samtools view -F 4 "$sam" | \
  awk '{
    seq=$10;
    len=length(seq);
    strand = (and($2,16)) ? "-" : "+";
    print seq "\t" len "\t" strand
  }'
done > LGTV_vsi_raw_all_updated.tsv # This is raw per-read information extracted from the SAM files. 

# Count identical reads
sort LGTV_vsi_raw_all_updated.tsv | uniq -c | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > LGTV_vsi_counts_all_updated.tsv
```

#### 2. Generate specific-stranded vsiRNA distribution Histogram

```sh
conda activate R_env
R
library(ggplot2)
library(readr)
library(dplyr)
library(scales)

# 2. Set working directory
setwd("/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific")

# 3. Read TSV file (no header)
vsi <- read_tsv(
  "LGTV_vsi_counts_all_updated.tsv",
  col_names = c("seq", "length", "strand", "count")
)
head(vsi)

# 4. Filter by length and extract 5' nucleotide (robust to lowercase)
vsi <- vsi %>%
  filter(length >= 18, length <= 30) %>%
  mutate(
    five_prime = toupper(substr(seq, 1, 1)),    # convert to uppercase
    five_prime = recode(five_prime, T = "U"),   # T -> U
    five_prime = factor(five_prime, levels = c("U", "A", "C", "G"))
  ) %>%
  filter(!is.na(five_prime))
head(vsi)

# 5. CPM normalization
total_mapped_reads <- sum(vsi$count)
> head(total_mapped_reads)
#[1] 274

vsi <- vsi %>%
  mutate(CPM = (count / total_mapped_reads) * 1e6)
head(vsi)

write_tsv(vsi, "LGTV_vsiRNA_CPM_per_sequence_updated.tsv")

# 6. Summarize CPM per length / strand / 5′ nucleotide
vsi_summary <- vsi %>%
  group_by(length, strand, five_prime) %>%
  summarise(CPM = sum(CPM), .groups = "drop")
head(vsi_summary)
write_tsv(vsi_summary, "LGTV_vsiRNA_CPM_summary_updated.tsv")

# 7. Create mirror values for minus strand
vsi_mirror <- vsi_summary %>%
  mutate(CPM_mirror = ifelse(strand == "-", -CPM, CPM))
head(vsi_mirror)
write_tsv(vsi_mirror, "LGTV_vsiRNA_CPM_mirror_updated.tsv")

# 8. Generate mirror plot and save as PDF
pdf("LGTV_vsiRNA_size_strand_5prime_CPM_mirror_updated.pdf", width = 6, height = 8)

ggplot(vsi_mirror, aes(x = length, y = CPM_mirror, fill = five_prime)) +
  geom_col(width = 0.8) +
  scale_x_continuous(breaks = 18:30) +
  scale_y_continuous(
    labels = function(x) comma(abs(x)),  # absolute + comma formatting
    name = "CPM"
  ) +
  scale_fill_manual(
    values = c("U" = "#E41A1C", "A" = "#4DAF4A",
               "C" = "#377EB8", "G" = "#984EA3")
  ) +
  labs(
    x = "vsiRNA length (nt)",
    fill = "5' nucleotide"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

dev.off()
```
- ##### Output Files
```sh
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific
```
