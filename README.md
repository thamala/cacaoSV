# cacaoSV
The repo contains the following software and scripts:
<br>
<br>
&nbsp;&nbsp;&nbsp;&nbsp;mumco_gaps.c: A program for filering SVs identified with MUM&Co based on their gap content<br>
&nbsp;&nbsp;&nbsp;&nbsp;mumco2vcf.c: A program for combining SVs indentified with MUM&Co into a joined VCF file<br>
&nbsp;&nbsp;&nbsp;&nbsp;SVstat.c: A program for estimating basic statistics (SFS, covariance matrix for PCA, *F*<sub>ST</sub>, *d*<sub>XY</sub>) for SVs<br>
&nbsp;&nbsp;&nbsp;&nbsp;SV2exp.c: A program for compiling RNA- and DNA-seq read counts for genes within 5 Kb of SVs<br>
&nbsp;&nbsp;&nbsp;&nbsp;SV2ase.c: A program for compiling RNA- and DNA-seq read counts for heterozygous SNPs found in heterozygous SVs<br>
&nbsp;&nbsp;&nbsp;&nbsp;SV2dxy.c: A program for estimating *d*<sub>XY</sub> and *F*<sub>ST</sub> between the SV arrangements<br>
&nbsp;&nbsp;&nbsp;&nbsp;SV2ld.c: A program for estimating LD within the major and minor SV arrangements<br>
&nbsp;&nbsp;&nbsp;&nbsp;SV2der.c A program for counting the number of derived alleles in the major and minor SV arrangements<br>
&nbsp;&nbsp;&nbsp;&nbsp;burn\_in.slim: A SLiM 3 script for conducting burn in for inv\_del.slim<br>
&nbsp;&nbsp;&nbsp;&nbsp;inv\_del.slim: A SLiM 3 script for simulating mutation accumulation inside inversions<br>
&nbsp;&nbsp;&nbsp;&nbsp;Exp\_lm.R: An R script for testing the association between SV genotypes and gene expression<br>
&nbsp;&nbsp;&nbsp;&nbsp;Admix\_lm.R: An R script for conducting an admixture-based genome scan<br>
