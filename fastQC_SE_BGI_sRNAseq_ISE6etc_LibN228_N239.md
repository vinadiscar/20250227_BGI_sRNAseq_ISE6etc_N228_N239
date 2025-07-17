# Quality Check of raw reads and clean reads using FastQC

- [FastQC](https://github.com/s-andrews/FastQC) is a program designed to spot potential problems in high througput sequencing datasets.
- In this analysis, we want to check the raw reads (beforming performing fastp) and check the clean reads (after performing fastp), then compare their quality.  By comparing the FastQC reports from both steps, we'll be able to assess how effective the cleaning process was.

##### Create directories to where fastQC data will be saved:

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC
mkdir raw # fastQC for raw reads (before fastp)
mkdir clean # fastQC for clean reads (after fastp)
```
## A. FastQC on Raw reads (before fastp)
- check the quality of raw reads (before cleaning with fastp)
  
### Write Script
```shell
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/raw
vi fastqc_raw.run
```
#### Script: "fastqc_raw.run"

```shell
#!/bin/bash
#SBATCH -t 4:00:00       # estimated runtime is less than 1 hour
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

        fastqc -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/raw \
         -t 22 ${QUERY}
done
```

### RUN fastqc_raw.run

```shell
sbatch fastqc_raw.run
# Submitted batch job 3510697 (Runtime: 25 mins. 45 secs.) 
```
### Unzip raw reads

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/raw
unzip '*.zip'
```


## B. FastQC on Clean reads (after fastp)

- Check the quality of clean reads (after cleaning with fastp)
  
### Write Script 

```shell
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/clean/
vi fastqc_clean.run
```
- #### Script: fastqc_clean.run

```shell
#!/bin/bash
#SBATCH -t 4:00:00       # estimated runtime is less than 1 hour
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
fastqc -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/clean -t 22 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/*.fq
```

### RUN fastqc_fastp_clean.run
```shell
sbatch fastqc_clean.run
# Submitted batch job 3510700 [Runtime: 03 mins., 04 secs.]
```

#### Unzip cleaned reads

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/clean
unzip '*.zip'
```
