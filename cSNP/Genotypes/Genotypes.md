---
title: WB MFM Genotypes per Gene
author: Deborah Velez-Irizarry
date: Fri Jan 11 09:54:15 EST 2019
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---


```r
### Description:  
```

Extract genotypes per gene from VCF files for Warmblood MFM project.  
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB    
  
File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;Genotypes.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB/Index_Chr_VariantCalling  
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB/Genotypes  

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;index_call_variants.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;index_vector.txt  
 
Render R Script  
  
> &nbsp;&nbsp;&nbsp;&nbsp;Index_Chr_VariantCalling.qsub  
  
***  
### R Environment    
Load required libraries


```r
library(vcfR)
```

```
## 
##    *****       ***   vcfR   ***       *****
##    This is vcfR 1.8.0 
##      browseVignettes('vcfR') # Documentation
##      citation('vcfR') # Citation
##    *****       *****      *****       *****
```

Clear Environment


```r
rm(list=ls())
```

Session Information


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_sandybridgep-r0.3.1.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] vcfR_1.8.0 knitr_1.20
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0        lattice_0.20-38   ape_5.1          
##  [4] viridisLite_0.3.0 permute_0.9-4     MASS_7.3-51.1    
##  [7] grid_3.5.1        nlme_3.1-137      magrittr_1.5     
## [10] evaluate_0.12     stringi_1.2.3     vegan_2.5-2      
## [13] Matrix_1.2-14     pinfsc50_1.1.0    tools_3.5.1      
## [16] stringr_1.3.1     parallel_3.5.1    compiler_3.5.1   
## [19] cluster_2.0.7-1   mgcv_1.8-24
```

### Load Data
VCF cSNP File


```r
vcf.Dir <- "/mnt/research/NMDL/MFM_WB/VariantCalling"
vcf.files <- list.files(vcf.Dir)[grep("vcf", list.files(vcf.Dir))]
genes <- unlist(lapply(strsplit(vcf.files, "_"), function(x) x[[1]][1]))
names(vcf.files) <- genes
```

Read VCF files into R


```r
vcf <- lapply(genes, function(x) read.vcfR(paste(vcf.Dir, vcf.files[x], sep="/"), verbose = FALSE))
names(vcf) <- genes
vcf
```

```
## $ACTA1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 17 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $BAG3
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 36 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $CRYAB
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 17 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $DES
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 57 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $DNAJB6
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 47 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $FHL1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 14 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $FLNC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 94 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $HSPB8
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 54 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $KY
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 292 variants
## Object size: 1.3 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $LMNA
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 17 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $MLIP
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 452 variants
## Object size: 1.2 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $MYOT
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 91 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PDLIM3
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 125 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PLEC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 112 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PLEC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 112 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PLEC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 112 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PLEC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 112 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PLEC
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 112 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $PYROXD1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 103 variants
## Object size: 1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $SQSTM1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 75 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $TIA1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 99 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
## 
## $TIAL1
## ***** Object of Class vcfR *****
## 16 samples
## 1 CHROMs
## 90 variants
## Object size: 1.1 Mb
## 0 percent missing data
## *****        *****         *****
```

Add ID column to VCF data by concating the chromosome and position of each cSNP


```r
vcf <- lapply(vcf, function(x) addID(x, sep = "_"))
head(vcf[[1]])
```

```
## [1] "***** Object of class 'vcfR' *****"
## [1] "***** Meta section *****"
## [1] "##fileformat=VCFv4.2"
## [1] "##FILTER=<ID=PASS,Description=\"All filters passed\">"
## [1] "##bcftoolsVersion=1.9-64-g28bcc56+htslib-1.9-52-g6e86e38"
## [1] "##bcftoolsCommand=mpileup -Ou -C50 -E -Q25 -a DV,AD,ADF,ADR,SP -f /m [Truncated]"
## [1] "##reference=file:///mnt/research/NMDL/MFM_WB/horse_r86_all_corrected.fa"
## [1] "##contig=<ID=10,length=83980604>"
## [1] "First 6 rows."
## [1] 
## [1] "***** Fixed section *****"
##      CHROM POS        ID           REF ALT QUAL  FILTER
## [1,] "1"   "68409137" "1_68409137" "T" "C" "199" "PASS"
## [2,] "1"   "68409294" "1_68409294" "T" "C" "85"  "PASS"
## [3,] "1"   "68409321" "1_68409321" "T" "C" "151" "PASS"
## [4,] "1"   "68409348" "1_68409348" "T" "C" "999" "PASS"
## [5,] "1"   "68409388" "1_68409388" "T" "C" "215" "PASS"
## [6,] "1"   "68409426" "1_68409426" "C" "T" "205" "PASS"
## [1] 
## [1] "***** Genotype section *****"
##      FORMAT                   12065.bam                         
## [1,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,135,255:0:0:29,0:16,0:45,0"
## [2,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,154,255:0:0:26,0:25,0:51,0"
## [3,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,151,255:0:0:23,0:27,0:50,0"
## [4,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,141,255:0:0:25,0:22,0:47,0"
## [5,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,135,255:0:0:19,0:26,0:45,0"
## [6,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,135,255:0:0:19,0:26,0:45,0"
##      12124.bam                          12130.bam                         
## [1,] "0/0:0,108,255:0:0:23,0:13,0:36,0" "0/0:0,131,255:1:0:31,0:21,1:52,1"
## [2,] "0/0:0,99,255:0:0:21,0:12,0:33,0"  "0/0:0,178,255:0:0:38,0:21,0:59,0"
## [3,] "0/0:0,111,255:0:0:26,0:11,0:37,0" "0/0:0,193,255:0:0:41,0:23,0:64,0"
## [4,] "1/1:255,99,0:33:0:0,20:0,13:0,33" "0/0:0,169,255:0:0:35,0:21,0:56,0"
## [5,] "0/0:0,114,255:0:0:18,0:20,0:38,0" "0/0:0,175,255:0:0:33,0:25,0:58,0"
## [6,] "0/0:0,96,255:0:0:14,0:18,0:32,0"  "0/0:0,160,255:0:0:28,0:25,0:53,0"
##      12212.bam                            
## [1,] "1/1:255,166,0:55:0:0,37:0,18:0,55"  
## [2,] "0/1:126,0,255:11:12:10,9:12,2:22,11"
## [3,] "0/1:191,0,255:18:7:11,13:11,5:22,18"
## [4,] "1/1:255,129,0:43:0:0,29:0,14:0,43"  
## [5,] "0/1:255,0,255:23:8:13,16:15,7:28,23"
## [6,] "0/0:0,160,255:0:0:26,0:27,0:53,0"   
##      12216.bam                         
## [1,] "0/0:0,102,255:0:0:26,0:8,0:34,0" 
## [2,] "0/0:0,87,255:0:0:22,0:7,0:29,0"  
## [3,] "0/0:0,93,255:0:0:22,0:9,0:31,0"  
## [4,] "1/1:255,102,0:34:0:0,25:0,9:0,34"
## [5,] "0/0:0,111,255:0:0:20,0:17,0:37,0"
## [6,] "0/0:0,108,255:0:0:17,0:19,0:36,0"
## [1] "First 6 columns only."
## [1] 
## [1] "Unique GT formats:"
## [1] "GT:PL:DV:SP:ADF:ADR:AD"
## [1]
```

### Review Genotypes
Extract genotype matrix per gene


```r
geno <- lapply(vcf, function(x) extract.gt(x, return.alleles=TRUE))

# Add animal names to column
anim <- unlist(lapply(strsplit(colnames(geno[[1]]), "[.]"), function(x) x[[1]][1]))
for (i in names(geno)){
    colnames(geno[[i]]) <- anim
}
```

Number of coding SNP per gene


```r
unlist(lapply(geno, nrow))
```

```
##   ACTA1    BAG3   CRYAB     DES  DNAJB6    FHL1    FLNC   HSPB8      KY 
##      17      36      17      57      47      14      94      54     292 
##    LMNA    MLIP    MYOT  PDLIM3    PLEC    PLEC    PLEC    PLEC    PLEC 
##      17     452      91     125     112     112     112     112     112 
## PYROXD1  SQSTM1    TIA1   TIAL1 
##     103      75      99      90
```

### Frequence of genotypes per gene


```r
freq <- lapply(geno, function(x) apply(x, 1, table))
```

### Save Genotype Information
Save called variant information to R data file


```r
save(vcf.files, vcf, geno, freq, file=paste(getwd(), "gene_variants.Rdata", sep="/"))
```

Save frequency per gene


```r
# Change genotype vector to data frame
freqData <- lapply(freq, function(x) lapply(x, function(y) 
    data.frame(rbind(genotypes=names(y), frequency=y))))

# Function to write genotype frequencies to file
add.snp <- function(gene, snp, gene.freq){
    write(snp, file=paste(gene, "freq.txt", sep="_"), append=TRUE)
    write.table(gene.freq, file=paste(gene, "freq.txt", sep="_"), 
        append=TRUE, col.names=FALSE, row.names=TRUE, quote=FALSE)
    write("", file=paste(gene, "freq.txt", sep="_"), append=TRUE)
}

# Create file per gene containing genotype frequency for all called coding SNP
for(i in names(freqData)){
    title <- paste("# Genotype frequencies per coding snp for", i, sep=" ")
    write(title, file=paste(i, "freq.txt", sep="_"))
    write("# SNP names consist of chromosome_position", file=paste(i, "freq.txt", sep="_"), append=TRUE)
    write("", file=paste(i, "freq.txt", sep="_"), append=TRUE)
    gene <- freqData[[i]]
    lapply(names(gene), function(x) add.snp(gene=i, snp=x, gene.freq=gene[[x]]))
}
```

Write genotype matrix per gene


```r
for (i in names(geno)){
    write.table(geno[[i]], file=paste(i, "genotypes.txt", sep="_"), 
        col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
}
```

### Run R Script


```r
htmlRunR
Genotypes.R nodes=1,cpus-per-task=1,time=03:00:00,mem=2G \
+WB MFM Genotypes per Gene
```

