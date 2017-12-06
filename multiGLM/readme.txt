The following data files and R scripts were used to produce results in the manuscript
"Multivariate Generalized Linear Model for Genetic Pleiotropy" by 
Daniel J. Schaid, Xingwei Tong, Anthony Batzler, Jason P. Sinnwell, 
Jiang Qing, Joanna M. Biernacka.


1) Scan SNPs with multivariate statistic t0 (testing null that no traits
   are associated with each SNP)

R Script:
    scan.multivariate.R


Data files: 
     dosage.fam:	PLINK fam file

     dosage.txt.gz:	SNP dosage file for 2001 SNPs on chromosome 15
     			note that only a prortion of chrom 15 is provided, because
			of the size of remaining chromosome files and	
			time to run scan.multivariate on all chromosomes.
			The SNP rs11635365 with the smallest p-value is this data.

     pleiotropy_phenotypes.txt:  phenotypes for 16 HAM-D questions.

     pleiotropy.15.results.txt:  output from scan.multivariate.R

2) Peform sequential sequential test of HAM-D response items with SNP rs11635365
   (Table 9 in manuscript).

R Script:
  sequential.R

Data files:
  sequential.dat.R:	text input to R as load("sequential.dat.R")

