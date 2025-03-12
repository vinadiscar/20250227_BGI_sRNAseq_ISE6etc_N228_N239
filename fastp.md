# Preprocess RNA seq data with fastp
[fastp](https://github.com/OpenGene/fastp) is an all-in-one tool for NGS read quality control.
- automatic adapter detection and removal
- N-containing read removal
- low-quality base filtering

## 1. Write Script [fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239.sbatch]

``` shell
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastp
vi fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239.sbatch
```
- #### Script: [fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239.sbatch]

``` shell
#!/bin/bash

#SBATCH -t 4:00:00       # Runtime - asking for 4  hours
#SBATCH -p cluster_short           # Partition (queue) - asking for large memory short queue
#SBATCH -J fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239             # Job name
#SBATCH -o %j.out            # Standard out
#SBATCH -e %j.err              # Standard error
#SBATCH -c 4
#SBATCH --mem 30G     # Memory needed per core
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

DIR1=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/
DIR2=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastp/

# Set up virtual environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh


# Preprocessing with fastp
conda activate

for LIB in N{228..239} ; do
        QUERY1_1=${DIR1}${LIB}/${LIB}"_L1_1.fq"
        FASTP_OUT1=${DIR2}"_L1_1.fq"

        fastp --in1 ${QUERY1_1} --out1 ${FASTP_OUT1} \
                  --qualified_quality_phred 30 --length_required 35 --overrepresentation_analysis --detect_adapter_for_pe \
                  -w 16 --html ${DIR2}.html --json ${DIR2}.json
done
```
## 2. Run [fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239.sbatch]

``` shell
sbatch fastp_20250227_BGI_sRNAseq_ISE6etc_N228-N239.sbatch
# Submitted job  [Runtime: 12 mins. 58 secs .]
```
### Output files: 
```sh
- /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess/fastp
```
