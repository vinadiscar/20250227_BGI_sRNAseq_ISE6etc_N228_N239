# Mapping Statistics of Total RNA seq ISE6 Mapped Against LGTV or HAZV genomes

### A. Mapped against LGTV genome
- ##### control samples
```sh
cd /work/ma-discar/20250128_BGI_hlo_ise_totalRNA/Analysis/STAR_Map_to_HAZV_LGTVgenomes_updated/STAR_Map_to_LGTVgenome
cat ISE6_NT_control_1_Log.final.out

Started job on |       Jan 30 19:24:57
                             Started mapping on |       Jan 30 19:24:57
                                    Finished on |       Jan 30 19:49:58
       Mapping speed, Million of reads per hour |       134.03

                          Number of input reads |       55883649
                      Average input read length |       286
                                    UNIQUE READS:
                   Uniquely mapped reads number |       13
                        Uniquely mapped reads % |       0.00%
                          Average mapped length |       264.92
                       Number of splices: Total |       0
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       0
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.41%
                         Deletion rate per base |       0.00%
                        Deletion average length |       0.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       0
             % of reads mapped to multiple loci |       0.00%
        Number of reads mapped to too many loci |       0
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       0.39%
                     % of reads unmapped: other |       99.61%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
- ##### LGTV-infected samples
```sh
cat ISE6_NT_LGTV_1_Log.final.out
Started job on |       Jan 30 18:59:47
                             Started mapping on |       Jan 30 18:59:47
                                    Finished on |       Jan 30 19:24:48
       Mapping speed, Million of reads per hour |       135.88

                          Number of input reads |       56652826
                      Average input read length |       285
                                    UNIQUE READS:
                   Uniquely mapped reads number |       506456
                        Uniquely mapped reads % |       0.89%
                          Average mapped length |       286.92
                       Number of splices: Total |       0
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       0
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.16%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       1.03
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       0
             % of reads mapped to multiple loci |       0.00%
        Number of reads mapped to too many loci |       169
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       0.34%
                     % of reads unmapped: other |       98.77%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
(base) ma-discar@cc21dev0:~/20250128_BGI_hlo_ise_totalRNA/Analysis/STAR_Map_to_HAZV_LGTVgenomes_updated/STAR_Map_to_LGTVgenome$ cat ISE6_NT_control_1_Log.final.out
                                 Started job on |       Jan 30 19:24:57
                             Started mapping on |       Jan 30 19:24:57
                                    Finished on |       Jan 30 19:49:58
       Mapping speed, Million of reads per hour |       134.03

                          Number of input reads |       55883649
                      Average input read length |       286
                                    UNIQUE READS:
                   Uniquely mapped reads number |       13
                        Uniquely mapped reads % |       0.00%
                          Average mapped length |       264.92
                       Number of splices: Total |       0
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       0
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.41%
                         Deletion rate per base |       0.00%
                        Deletion average length |       0.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       0
             % of reads mapped to multiple loci |       0.00%
        Number of reads mapped to too many loci |       0
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       0.39%
                     % of reads unmapped: other |       99.61%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```


### B. Mapped against HAZV genome 

- ##### control samples
```sh
cd /work/ma-discar/20250128_BGI_hlo_ise_totalRNA/Analysis/STAR_Map_to_HAZV_LGTVgenomes_updated/STAR_Map_to_HAZVgenome
cat ISE6_NT_control_2_Log.final.out
 Started job on |       Jan 30 19:37:04
                             Started mapping on |       Jan 30 19:37:04
                                    Finished on |       Jan 30 19:52:54
       Mapping speed, Million of reads per hour |       96.90

                          Number of input reads |       25569917
                      Average input read length |       286
                                    UNIQUE READS:
                   Uniquely mapped reads number |       24
                        Uniquely mapped reads % |       0.00%
                          Average mapped length |       196.29
                       Number of splices: Total |       0
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       0
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.83%
                         Deletion rate per base |       0.00%
                        Deletion average length |       0.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       0
             % of reads mapped to multiple loci |       0.00%
        Number of reads mapped to too many loci |       0
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       3.56%
                     % of reads unmapped: other |       96.44%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
- ##### HAZV-infected samples 
```sh
cat ISE6_NT_HAZV_1_Log.final.out
Started job on |       Jan 30 18:58:38
                             Started mapping on |       Jan 30 18:58:38
                                    Finished on |       Jan 30 19:36:54
       Mapping speed, Million of reads per hour |       95.11

                          Number of input reads |       60658315
                      Average input read length |       283
                                    UNIQUE READS:
                   Uniquely mapped reads number |       4826614
                        Uniquely mapped reads % |       7.96%
                          Average mapped length |       286.14
                       Number of splices: Total |       0
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       0
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.18%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       1.09
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       0
             % of reads mapped to multiple loci |       0.00%
        Number of reads mapped to too many loci |       7299
             % of reads mapped to too many loci |       0.01%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       3.89%
                     % of reads unmapped: other |       88.14%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
