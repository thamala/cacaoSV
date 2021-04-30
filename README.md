# cacaoSV
The repo contains the following software and scripts:
<br>
<br>
mumco_gaps.c: A program for filtering SVs identified with MUM&Co based on their gap content<br>
mumco2vcf.c: A program for combining SVs indentified with MUM&Co into a joined VCF file<br>
SVstat.c: A program for estimating basic statistics (SFS, covariance matrix for PCA, *F*<sub>ST</sub>, *d*<sub>XY</sub>) for SVs<br>
SV2exp.c: A program for compiling RNA- and DNA-seq read counts for genes within 5 Kb of SVs<br>
SV2ase.c: A program for compiling RNA- and DNA-seq read counts for heterozygous SNPs in heterozygous SVs<br>
SV2dxy.c: A program for estimating *d*<sub>XY</sub> and *F*<sub>ST</sub> between the SV arrangements<br>
SV2ld.c: A program for estimating LD within the major and minor SV arrangements<br>
SV2der.c A program for counting the number of derived alleles in the major and minor SV arrangements<br>
SV2sift: A program for compiling SIFT scores for SNPs in the major and minor SV arrangements<br>
burn\_in.slim: A SLiM 3 script for conducting burn in for inv\_del.slim<br>
inv\_del.slim: A SLiM 3 script for simulating mutation accumulation inside inversions<br>
Exp\_lm.R: An R script for testing the association between SV genotypes and gene expression<br>
Admix\_lm.R: An R script for conducting an admixture-based genome scan<br>
