## script to run t0 test (multivariate test of null hypothesis that no traits
## are associated with a SNP) for all SNP's on chromosome 15

## Set input parameters
chr=15
subset=""
pheno<-"pleiotropy"

### read in SNP data
dosage<-read.table(file="dosage.txt.gz",header=T,row.names=1,check.names=F)
alleles<-dosage[,1:2]
dosage<-dosage[,-c(1:2)]
## subset dosage if subset provided
fam<-read.table(file="dosage.fam",header=F)
if(subset!=""){
  if(subset=="male"){
    fam<-fam[fam[,5]==1,]
  } else if(subset=="female"){
    fam<-fam[fam[,5]==2,]
  }
} 
dosage<-dosage[,as.character(fam[,2])]
header<-colnames(dosage)

## Get phenotype and dosage on same subjects
b.hamd.q<-read.table(file="pleiotropy_phenotypes.txt",header=T)
keep<-b.hamd.q$IID[b.hamd.q$IID %in% header]
b.hamd.q<-b.hamd.q[b.hamd.q$IID %in% keep,]
dosage<-dosage[,as.character(keep)]

## Now that we've subset the dosage, recalculate MAF, remove
## sb.hamd.qs where MAF < .01 and flip alleles and dosage so based on new MAF
MAF<-round(apply(dosage,1,function(x) sum(x,na.rm=T)/(2*sum(!is.na(x)))),3)

## remove any where MAF < .01 (note: could be coded by major or minor allele)
remove<-MAF < .01 | MAF > .99
dosage<-dosage[!remove,]
MAF<-MAF[!remove]
alleles<-alleles[!remove,]

## Flip dosage where MAF > .5
flip<-MAF > .5
dosage[flip,]<-2-dosage[flip,]
MAF[flip]<-1-MAF[flip]
ma<-ifelse(flip,alleles[,1],alleles[,2])
ca<-ifelse(flip,alleles[,2],alleles[,1])
rm(alleles)

dosage.other<-cbind(ma,ca,MAF)
  

## Note that plieo package (from CRAN) should be installed
library(pleio)


## make matrix of only responses to questions on HAMD  ###
y.mat<-as.matrix(b.hamd.q[,-c(1:2)])

## need to define the type of family for traits:
glm.family <- rep("gaussian", ncol(y.mat))

pval0 <- rep(NA, nrow(dosage))

for(i in 1:nrow(dosage)){
  snp <-dosage[i,]
  fit <-  pleio.glm.fit(y.mat, snp, glm.family)
  pval0[i]<- pleio.glm.test(fit, count.nonzero.coef = 0)$pval
}

tbl <- cbind(format(pval0,digits=4,scientific=T),dosage.other)
write.table(tbl, file=paste(pheno,subset,".",chr,".results.txt",sep=""),row.names=T,col.names=F,quote=F)
