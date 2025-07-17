# Mapping of cleaned 20250227_BGI_sRNAseq_ISE6etc RNA seq data (H longicornis) to HaeL2018 Genome
- Map cleaned 20250227_BGI_sRNAseq_ISE6etc RNA seq data (Okayama samples, N228 to N239) to HaeL2018 genome using STAR
-  HaeL2018 genome was downloaded from [Vectorbase](https://vectorbase.org/vectorbase/app)
-  Process: I. Indexing, II. Mapping
- Tailored to Single End reads and smallRNAseq Data


## I. Indexing

```sh
conda activate star-env
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR
mkdir index # create folder where index results will be found
```

- ### Write Script: [STAR_GenomeIndex_HaeL2018.run]

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR/index
vi STAR_GenomeIndex_HaeL2018.run
```
##### Script: [STAR_GenomeIndex_HaeL2018.run]

```sh
#!/bin/bash

#SBATCH -p cluster_short           # partition name
#SBATCH -t 4:00:00              # hours:minutes runlimit after which job will be killed
#SBATCH -n 30           # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 80G
#SBATCH --job-name STAR_index           # Job name
#SBATCH -o %j.out                       # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH--mail-type=ALL

# Set up virtual environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh

# StarIndex
conda activate star-env

STAR --runThreadN 30 \
--runMode genomeGenerate \
--genomeDir /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR/index \
--genomeFastaFiles /work/ma-discar/H_longicornis/GenomeInfo_HaeL2018/VectorBase-66_HlongicornisHaeL2018_Genome.fasta \
--sjdbGTFfile /work/ma-discar/H_longicornis/GenomeInfo_HaeL2018/VectorBase-66_HlongicornisHaeL2018.gtf \
--sjdbOverhang 99

```

- ### Run: [STAR_GenomeIndex_HaeL2018.run]

```sh
sbatch STAR_GenomeIndex_HaeL2018.run
# Submitted batch job 3519691 (Runtime: 2 hrs., 33 mins.)
```

## II. Mapping of cleaned RNA Seq data to HaeL208 genome

- ### Write Script

```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR
vi STAR_Map_HaeL2018.run
```

##### Script [STAR_Map_HaeL2018.run]

```sh
#!/bin/bash

#SBATCH -p gpu_short                 # partition name
#SBATCH -t 4:00:00                       # hours:minutes runlimit after which job will be killed
#SBATCH -n 4                             # number of cores requested
#SBATCH --mem 0                          # grants access to all of the memory on each node
#SBATCH --job-name STAR_map             # Job name
#SBATCH -o %j.out                        # File to which standard out will be written
#SBATCH -e %j.err                        # File to which standard err will be written
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

# Activate conda environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh
conda activate star-env

# Define directories
DIR1=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp
DIR2=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR
GENOME_INDEX=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR/index

# STAR Mapping
for i in {228..239}; do
    LIB="N${i}_L1"
    FASTP_OUT1="${DIR1}/${LIB}_1.fq"

    STAR --runThreadN 4 \
         --genomeDir ${GENOME_INDEX} \
         --readFilesIn ${FASTP_OUT1} \
         --outFileNamePrefix ${DIR2}/HaeL_${LIB} \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 20 \
         --outFilterMismatchNmax 2 \
         --outFilterMismatchNoverReadLmax 0.05 \
         --alignIntronMax 1 \
         --outFilterScoreMinOverLread 0.66 \
         --outFilterMatchNminOverLread 0.66
done

# Output file: ***Aligned.sortedByCoord.out.bam
```

- ### Run

```sh
sbatch STAR_Map_HaeL2018.run
#Submitted batch job 3529874  (Runtime:2hrs:39secs )
```

## III. Check Mapping Statistics

* Checking alignment statistics in ".Log.final.out" files using multiQC:

- #### Write Script [multiqc_StarHaeL2018.sh]
```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR
mkdir multiQc_mapStat 
vi multiqc_StarHaeL2018.sh
```
- multiqc_StarHaeL2018.sh:

```sh
#!/bin/bash

# Define the directory containing STAR Log.final.out files
log_dir="/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR"  # Change the directory where Log.final.out files are stored

# Run MultiQC on the directory containing Log.final.out files
echo "Running MultiQC on $log_dir..."

multiqc $(find $log_dir -type f -name "*Log.final.out") -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/Mapping_to_HaeL2018genome_STAR/multiQc_mapStat/multiqc_StarHaeL2018.html

# Output report location
echo "MultiQC report generated in multiqc_report.html"
```

- #### Run the script multiqc_StarHaeL2018.sh
```sh
conda activate multiqc_env
chmod +x multiqc_StarHaeL2018.sh # make the script executable 
./multiqc_StarHaeL2018.sh # run the script
```
