# Map small RNAseq data Hlo to HAZV genome using Bowtie 1
- Bowtie 1 (rather than Bowtie 2) is recommended for small RNA-seq data because Bowtie 1 is optimized for short reads (e.g., 18–30 nt), typical of small RNAs.

- ##### Install Bowtie 1 using Conda
> conda create -n bowtie1_env -c bioconda bowtie=1.2.3 -y
> conda activate bowtie1_env

## 1. Build the Bowtie Index
```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1
mkdir bowtie_index
```
- ### Write script: [bowtie1_index_sRNAseq_HAZV.run]
```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/bowtie_index
vi bowtie1_index_sRNAseq_HAZV.run
```
- ###### bowtie1_index_sRNAseq_HAZV.run:
```sh
#!/bin/bash

#SBATCH -p cluster_short           # Partition name
#SBATCH -t 4:00:00                 # Walltime limit
#SBATCH -n 12                      # Number of CPU cores
#SBATCH --mem=80G                  # Total memory
#SBATCH --job-name=bowtie1_HAZVgenome_index
#SBATCH -o %j.out                  # STDOUT log (job ID as filename)
#SBATCH -e %j.err                  # STDERR log (job ID as filename)
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

# Load conda environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh
conda activate bowtie1_env

# Change to the working directory where output should go 
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/bowtie_index

# Run Bowtie 1 index builder
bowtie-build /work/ma-discar/Viral_genome_info_NCBI/HAZV_ncbi_dataset/GCA_002831085.1_ASM283108v1_genomic.fna hazv_index

```
- ### Run script: [bowtie1_index_sRNAseq_HAZV.run]
```sh
sbatch bowtie1_index_sRNAseq_HAZV.run
# submitted job ID 3510704  (Runtime: 5 seconds)
```
## 2. Map HAZV-infected samples against HAZV genome

- ### Write Script [bowtie1_map_sRNAseq_HAZV.run]

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result
vi bowtie1_map_SE_HlosRNAseq_HAZV.run
```
- ###### bowtie1_map_SE_HlosRNAseq_HAZV.run
```sh
#!/bin/bash

#SBATCH -p cluster_short
#SBATCH -t 4:00:00
#SBATCH -n 4
#SBATCH --mem=40G
#SBATCH --job-name=Bowtie1_smallRNAseq_HAZV
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

# Activate conda environment with Bowtie1
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh
conda activate bowtie1_env

# Define directories
DIR1=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/
DIR2=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/mapping_result/
GENOME_INDEX=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HAZVgenome_bowtie1/bowtie_index/hazv_index

# Create output directory if it doesn't exist
mkdir -p ${DIR2}/logs

# Process all samples N228–N239
for i in {228..239}; do
    LIB="N${i}_L1"
    READ1="${DIR1}/${LIB}_1.fq"
    OUT_SAM="${DIR2}/HAZVgenome_${LIB}.sam"
    LOG_FILE="${DIR2}/logs/${LIB}_bowtie.log"

    echo "Mapping ${LIB} to HAZV genome..."
    bowtie -v 1 -k 10 --best --strata -p 4 \
        ${GENOME_INDEX} \
        ${READ1} \
        ${OUT_SAM} 2> ${LOG_FILE}

    echo "Log saved to ${LOG_FILE}"
done

```
## Run Script
```sh
sbatch bowtie1_map_SE_HlosRNAseq_HAZV.run
# Submitted batch job ID 3510712 (Runtime:04 mins., 37 secs.)
```
