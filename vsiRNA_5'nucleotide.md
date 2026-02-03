# Strand-specific vsiRNA size distribution: 5'nucleotide Analysis

Mapping to viral genome (Bowtie1)
  ↓
SAM file
  ↓
Count vsiRNAs (awk/samtools)
  ↓
vsi_counts_all.tsv   ← raw counts
  ↓
CPM normalization
  ↓
Plots/figures


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
## II HAZV-Infected (Update): Each Day/Time Point
#### 1. Generate Input Frame (Extract Raw Reads)
```sh
conda activate samtools_env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/vsiRNA_Profiles/vsiRNA_profile_Strandspecific/vsiRNA_5nucleotide_TimePoints

# Declare SAM files and corresponding day
declare -A sam_day
sam_day=(
  [N230]="/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N230_L1.sam"
  [N233]="/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N233_L1.sam"
  [N236]="/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N236_L1.sam"
  [N239]="/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result_update/HAZVgenome_N239_L1.sam"
)

for day in "${!sam_day[@]}"; do
    sam="${sam_day[$day]}"
    
    echo "Processing $day from $sam ..."

    # Extract mapped reads
    samtools view -F 4 "$sam" | \
    awk '{
        ref=$3; seq=$10; len=length(seq);      
        strand = (and($2,16)) ? "-" : "+";
        print ref "\t" seq "\t" len "\t" strand
    }' > "HAZV_vsi_raw_${day}.tsv"

    # Collapse identical sequences
    LC_ALL=C sort "HAZV_vsi_raw_${day}.tsv" | uniq -c | \
    awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1}' \
    > "HAZV_vsi_counts_${day}.tsv"
done
```
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
