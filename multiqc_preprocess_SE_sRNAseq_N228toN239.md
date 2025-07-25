# Aggregate Preprocessing Results with MultiQC
- This method summarizes preprocessing results using multiQC
## Multiqc
- [Multiqc](https://multiqc.info/) is a modular tool to aggregate results from bioinformatics analyses across many samples into a single report.
- It searches a given directory for analysis logs and compiles a HTML report.
- It's a general use tool, perfect for summarising the output from numerous bioinformatics tools

### Activate multiQC environment in Conda
```shell
cd /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/
conda activate multiqc_env
```

## Run MultiQC

### A. Integrate fastqc (raw) reports
```shell
multiqc /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/raw -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/multiqc_report_raw -n multiqc_report_raw.html
```
### B. Integrate fastqc (clean) reports
```shell
multiqc /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/fastQC/clean -o /work/ma-discar/20250227_BGI_sRNAseq_ISE6etc_N228_N239/preprocess_update/multiqc_report_clean -n multiqc_report_clean.html
```

