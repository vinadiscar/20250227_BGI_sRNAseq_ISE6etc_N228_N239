# Preprocess small RNAseq data with fastp
- preprocess small RNA seq data (BGI_sRNAseq_ISE6etc, Libraries N228 to N239)
- Single End

## 1. Write Script
```sh
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp

vi fastp_SE_Hlo_sRNAseq.sbatch
```
- fastp_SE_Hlo_sRNAseq.sbatch

```sh
#!/bin/bash

#SBATCH -t 4:00:00       # Runtime - asking for 4  hours
#SBATCH -p cluster_short           # Partition (queue) - asking for large memory short queue
#SBATCH -J fastp_SE_HlosRNAseq_ISE6etc_N228-N239             # Job name
#SBATCH -o %j.out            # Standard out
#SBATCH -e %j.err              # Standard error
#SBATCH -c 4
#SBATCH --mem 30G     # Memory needed per core
#SBATCH --mail-user=discar.ma_divina_kristi.di5@naist.ac.jp
#SBATCH --mail-type=ALL

DIR1=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/
DIR2=/work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/

# Set up virtual environment
source /work/ma-discar/anaconda3/etc/profile.d/conda.sh


# Preprocessing with fastp
conda activate

for LIB in N{228..239} ; do
        QUERY1_1=${DIR1}${LIB}/${LIB}"_L1_1.fq"
        FASTP_OUT1=${DIR2}${LIB}"_L1_1.fq"

        fastp \
  --in1 ${QUERY1_1} --out1 ${FASTP_OUT1} \
  --adapter_sequence AGATCGGAAGAGCACACGTCT \
  --qualified_quality_phred 30 \
  --length_required 18 \
  --overrepresentation_analysis \
  -w 4 \
  --html ${DIR2}${LIB}_L1_1.html \
  --json ${DIR2}${LIB}_L1_1.json

done
```
## 2. Run script [fastp_SE_Hlo_sRNAseq.sbatch]
```sh
sbatch fastp_SE_Hlo_sRNAseq.sbatch
# Submitted batch job 3510686 (Runtime: 11 mins., 17 secs.)
```

-----
## Log:
```sh
Read1 before filtering:
total reads: 44406887
total bases: 2220344350
Q20 bases: 2180063224(98.1858%)
Q30 bases: 2088804121(94.0757%)

Read1 after filtering:
total reads: 42821747
total bases: 1238687627
Q20 bases: 1223408960(98.7665%)
Q30 bases: 1187089071(95.8344%)

Filtering result:
reads passed filter: 42821747
reads failed due to low quality: 932816
reads failed due to too many N: 1
reads failed due to too short: 652323
reads with adapter trimmed: 38955927
bases trimmed due to adapters: 944055108

Duplication rate (may be overestimated since this is SE data): 75.1304%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N228_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N228_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N228/N228_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N228_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N228_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N228_L1_1.json
fastp v0.23.2, time used: 48 seconds
Read1 before filtering:
total reads: 53797304
total bases: 2689865200
Q20 bases: 2642484317(98.2385%)
Q30 bases: 2537209490(94.3248%)

Read1 after filtering:
total reads: 51502614
total bases: 1388791669
Q20 bases: 1372606272(98.8346%)
Q30 bases: 1334507183(96.0912%)

Filtering result:
reads passed filter: 51502614
reads failed due to low quality: 1203981
reads failed due to too many N: 3
reads failed due to too short: 1090706
reads with adapter trimmed: 49685167
bases trimmed due to adapters: 1250843678

Duplication rate (may be overestimated since this is SE data): 80.5207%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N229_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N229_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N229/N229_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N229_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N229_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N229_L1_1.json
fastp v0.23.2, time used: 57 seconds
Read1 before filtering:
total reads: 59624032
total bases: 2981201600
Q20 bases: 2922075527(98.0167%)
Q30 bases: 2790597303(93.6065%)

Read1 after filtering:
total reads: 56956958
total bases: 1526232365
Q20 bases: 1506689452(98.7195%)
Q30 bases: 1460777147(95.7113%)

Filtering result:
reads passed filter: 56956958
reads failed due to low quality: 1455251
reads failed due to too many N: 2
reads failed due to too short: 1211821
reads with adapter trimmed: 54530992
bases trimmed due to adapters: 1395689069

Duplication rate (may be overestimated since this is SE data): 79.8532%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N230_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N230_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N230/N230_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N230_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N230_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N230_L1_1.json
fastp v0.23.2, time used: 62 seconds
Read1 before filtering:
total reads: 62844902
total bases: 3142245100
Q20 bases: 3096391419(98.5407%)
Q30 bases: 2993966954(95.2811%)

Read1 after filtering:
total reads: 59081380
total bases: 1575386553
Q20 bases: 1559224972(98.9741%)
Q30 bases: 1520656681(96.5259%)

Filtering result:
reads passed filter: 59081380
reads failed due to low quality: 1365765
reads failed due to too many N: 9
reads failed due to too short: 2397748
reads with adapter trimmed: 58135459
bases trimmed due to adapters: 1494418524

Duplication rate (may be overestimated since this is SE data): 79.2045%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N231_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N231_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N231/N231_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N231_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N231_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N231_L1_1.json
fastp v0.23.2, time used: 65 seconds
Read1 before filtering:
total reads: 62752988
total bases: 3137649400
Q20 bases: 3087245423(98.3936%)
Q30 bases: 2975210751(94.8229%)

Read1 after filtering:
total reads: 56595501
total bases: 1568636405
Q20 bases: 1550914383(98.8702%)
Q30 bases: 1509291864(96.2168%)

Filtering result:
reads passed filter: 56595501
reads failed due to low quality: 1288524
reads failed due to too many N: 5
reads failed due to too short: 4868958
reads with adapter trimmed: 57811678
bases trimmed due to adapters: 1461715662

Duplication rate (may be overestimated since this is SE data): 80.9046%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N232_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N232_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N232/N232_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N232_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N232_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N232_L1_1.json
fastp v0.23.2, time used: 67 seconds
Read1 before filtering:
total reads: 55102634
total bases: 2755131700
Q20 bases: 2705549167(98.2004%)
Q30 bases: 2593993030(94.1513%)

Read1 after filtering:
total reads: 49642452
total bases: 1349921793
Q20 bases: 1333762914(98.803%)
Q30 bases: 1295401509(95.9612%)

Filtering result:
reads passed filter: 49642452
reads failed due to low quality: 1192009
reads failed due to too many N: 2
reads failed due to too short: 4268171
reads with adapter trimmed: 50928809
bases trimmed due to adapters: 1307628910

Duplication rate (may be overestimated since this is SE data): 76.8579%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N233_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N233_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N233/N233_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N233_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N233_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N233_L1_1.json
fastp v0.23.2, time used: 59 seconds
Read1 before filtering:
total reads: 34588629
total bases: 1729431450
Q20 bases: 1674318721(96.8132%)
Q30 bases: 1568289961(90.6824%)

Read1 after filtering:
total reads: 32408828
total bases: 871596756
Q20 bases: 851599476(97.7057%)
Q30 bases: 812397422(93.2079%)

Filtering result:
reads passed filter: 32408828
reads failed due to low quality: 1706736
reads failed due to too many N: 4
reads failed due to too short: 473061
reads with adapter trimmed: 31436644
bases trimmed due to adapters: 802867126

Duplication rate (may be overestimated since this is SE data): 74.4253%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N234_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N234_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N234/N234_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N234_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N234_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N234_L1_1.json
fastp v0.23.2, time used: 38 seconds
Read1 before filtering:
total reads: 41557149
total bases: 2077857450
Q20 bases: 2013919976(96.9229%)
Q30 bases: 1888874123(90.9049%)

Read1 after filtering:
total reads: 39009253
total bases: 1134651576
Q20 bases: 1108351302(97.6821%)
Q30 bases: 1056423077(93.1055%)

Filtering result:
reads passed filter: 39009253
reads failed due to low quality: 1806873
reads failed due to too many N: 4
reads failed due to too short: 741019
reads with adapter trimmed: 36134723
bases trimmed due to adapters: 878177084

Duplication rate (may be overestimated since this is SE data): 70.039%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N235_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N235_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N235/N235_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N235_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N235_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N235_L1_1.json
fastp v0.23.2, time used: 48 seconds
Read1 before filtering:
total reads: 44742509
total bases: 2237125450
Q20 bases: 2182222553(97.5458%)
Q30 bases: 2076587380(92.8239%)

Read1 after filtering:
total reads: 41075993
total bases: 1129777801
Q20 bases: 1108103235(98.0815%)
Q30 bases: 1065098553(94.275%)

Filtering result:
reads passed filter: 41075993
reads failed due to low quality: 1874016
reads failed due to too many N: 6
reads failed due to too short: 1792494
reads with adapter trimmed: 40948151
bases trimmed due to adapters: 1027881803

Duplication rate (may be overestimated since this is SE data): 71.3445%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N236_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N236_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N236/N236_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N236_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N236_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N236_L1_1.json
fastp v0.23.2, time used: 50 seconds
Read1 before filtering:
total reads: 58300016
total bases: 2915000800
Q20 bases: 2832903576(97.1836%)
Q30 bases: 2673490772(91.7149%)

Read1 after filtering:
total reads: 54291359
total bases: 1486782548
Q20 bases: 1455367105(97.887%)
Q30 bases: 1392990058(93.6916%)

Filtering result:
reads passed filter: 54291359
reads failed due to low quality: 2696598
reads failed due to too many N: 6
reads failed due to too short: 1312053
reads with adapter trimmed: 52883073
bases trimmed due to adapters: 1333408780

Duplication rate (may be overestimated since this is SE data): 71.9362%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N237_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N237_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N237/N237_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N237_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N237_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N237_L1_1.json
fastp v0.23.2, time used: 64 seconds
Read1 before filtering:
total reads: 46309472
total bases: 2315473600
Q20 bases: 2261569890(97.672%)
Q30 bases: 2157935274(93.1963%)

Read1 after filtering:
total reads: 43529936
total bases: 1208662258
Q20 bases: 1186128451(98.1356%)
Q30 bases: 1141116993(94.4116%)

Filtering result:
reads passed filter: 43529936
reads failed due to low quality: 1944090
reads failed due to too many N: 3
reads failed due to too short: 835443
reads with adapter trimmed: 42007115
bases trimmed due to adapters: 1039877175

Duplication rate (may be overestimated since this is SE data): 72.5593%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N238_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N238_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N238/N238_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N238_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N238_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N238_L1_1.json
fastp v0.23.2, time used: 50 seconds
Read1 before filtering:
total reads: 60696969
total bases: 3034848450
Q20 bases: 2931226596(96.5856%)
Q30 bases: 2733095541(90.0571%)

Read1 after filtering:
total reads: 56875459
total bases: 1569584041
Q20 bases: 1530648412(97.5194%)
Q30 bases: 1455310329(92.7195%)

Filtering result:
reads passed filter: 56875459
reads failed due to low quality: 3049089
reads failed due to too many N: 6
reads failed due to too short: 772415
reads with adapter trimmed: 54798809
bases trimmed due to adapters: 1364804109

Duplication rate (may be overestimated since this is SE data): 70.4641%

JSON report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N239_L1_1.json
HTML report: /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N239_L1_1.html

fastp --in1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/N239/N239_L1_1.fq --out1 /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N239_L1_1.fq --adapter_sequence AGATCGGAAGAGCACACGTCT --qualified_quality_phred 30 --length_required 18 --overrepresentation_analysis -w 4 --html /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N239_L1_1.html --json /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastp/N239_L1_1.json
fastp v0.23.2, time used: 67 seconds
```
