# Quality Check of raw reads and clean reads using FastQC

- [FastQC](https://github.com/s-andrews/FastQC) is a program designed to spot potential problems in high througput sequencing datasets.
- In this analysis, we want to check the raw reads (beforming performing fastp) and check the clean reads (after performing fastp), then compare their quality.  By comparing the FastQC reports from both steps, we'll be able to assess how effective the cleaning process was in improving our data's quality.

##### Create directories to where fastQC data will be saved:

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc
mkdir raw # fastQC for raw reads (before fastp)
mkdir clean # fastQC for clean reads (after fastp)
```
## A. FastQC on Raw reads (before fastp)
- check the quality of raw reads (before fastp)
  
### Write Script for FastQC on raw reads(before fastp)
```shell
/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/raw
vi fastqc_raw.run
```
#### Script: "fastqc_raw.run"

```shell
#!/bin/bash
#SBATCH -t 4:00:00       # Maximum runtime; estimated runtime is less than 1 hour
#SBATCH -p cluster_short # Partition (queue) - asking for large memory short queue
#SBATCH -J fastqc_raw    # Job name
#SBATCH -o %j.out            # Standard out
#SBATCH -e %j.err              # Standard error
#SBATCH -c 22          # number of cores
#SBATCH --mem 30G     # Memory needed per core
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

## Set up virtual environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh

## ActivateConda Environment for FASTQC
conda activate myenv

DIR=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/

for LIB in N{228..239} ; do
        QUERY=${DIR}${LIB}/${LIB}"_L1_1.fq"

        fastqc -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/raw \
         -t 22 ${QUERY}
done

```

### RUN fastqc_raw.run

```shell
sbatch fastqc_raw.run
# Submitted batch job 3463674(Runtime: 25 mins., 21 secs.) 
```
### Unzip raw reads

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/raw
unzip '*.zip'
```


## B. FastQC on Clean reads (after fastp)

- Check the quality of clean reads (after fastp)
  
### Write Script for FastQC (after fastp)

```shell
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/clean/
vi fastqc_fastp_clean.run
```
#### fastqc_fastp_clean.run

```shell
#!/bin/bash
#SBATCH -t 1:00:00       # Maximum runtime; estimated runtime is less than 1 minute
#SBATCH -p cluster_short # Partition (queue) - asking for large memory short queue
#SBATCH -J fastqc_clean    # Job name
#SBATCH -o %j.out            # Standard out
#SBATCH -e %j.err              # Standard error
#SBATCH -c 22          # number of cores
#SBATCH --mem 30G     # Memory needed per core
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

## Set up virtual environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh

## ActivateConda Environment for FASTQC
conda activate myenv

## Run FASTQC
fastqc -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/clean -t 22 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastp/*.fq
```

### RUN fastqc_fastp_clean.run
```shell
sbatch fastqc_fastp_clean.run
# Submitted batch job 3463677 [Runtime: 34 secs.]
```

#### Unzip cleaned reads

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastqc/clean
unzip '*.zip'
```
