#Load libraries
library(qqman)
library(data.table)
library(dplyr)
# Attaching package: 'dplyr'
## The following objects are masked from 'package:data.table':
## 
##     between, first, last
## The following objects are masked from 'package:data.table':
## 
##     between, first, last
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union

#paste function to concatenate filenames/paths
"%&%" = function(a,b) paste(a,b,sep="")
#paths to gwas sumstat files
data.dir <- "/home/vir/TypeII/"
T2 <- data.dir %&% "phecode-250.2-both_sexes.tsv.bgz"
## command to read in a .bgz zipped file
gwas <- fread(cmd="zcat " %&% T2)
#replace X chr with 23 to plot
gwas <- mutate(gwas,chr=as.numeric(ifelse(chr=="X",23,chr)))
gwas <- mutate(gwas, p=10^(-1 * neglog10_pval_CSA))
#newpvals<-select(gwas, chr, pos, p)%>%arrange(p)
#newpvals[1:40,3]<-1e-30
#newgwas<-left_join(gwas,newpvals,by=(c("chr","pos")))
newgwas <- mutate(gwas, p=base::ifelse(p==0, 1e-30, p))
# Check if else statement^^
# Different types of ifelse --> library name :: function
#filter meta results to SNPs with P<0.01 for fewer dots to plot
gwasREAL <- dplyr::filter(newgwas,p < 0.01)
gwasREAL <- dplyr::filter(newgwas, af_controls_CSA > 0.01 & af_controls_CSA < 0.99)



manhattan(gwasREAL,chr="chr",bp="pos",p="p",snp="pos",suggestiveline = FALSE,
          main="Pan-UK Biobank Type II Diabetes CSA (n=8,845)")

## BELOW: GDM Cases to Controls Per Pop Plotted
## CSA --> 75 cases to 444 controls
## EUR -->  811 cases to 6432 controls
## META --> 886 cases to 6876 controls

## BELOW: Type II Cases to Controls Per Pop Plotted
## CSA --> 1664 cases to 7181 controls
## EUR --> 22768 cases to 396181 controls
## META --> 255550 cases to 413134 controls
## AFR --> 796 cases to 5793 controls
## EAS --> 152 cases to 2555 controls
## MID --> 170 cases to 1424 controls
