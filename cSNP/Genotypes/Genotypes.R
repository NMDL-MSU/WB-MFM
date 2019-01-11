### Description:  
#' Extract genotypes per gene from VCF files for Warmblood MFM project.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB    
#'   
#' File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;Genotypes.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB/Index_Chr_VariantCalling  
#'   
#' **Output files:**  
#'   
#' Directory:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/MFM_WB/Genotypes  
#' 
#' Files:
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;index_call_variants.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;index_vector.txt  
#'  
#' Render R Script  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;Index_Chr_VariantCalling.qsub  
#'   
#' ***  

#' ### R Environment    

#' Load required libraries
library(vcfR)

#' Clear Environment
rm(list=ls())

#' Session Information
sessionInfo()


#' ### Load Data

#' VCF cSNP File
vcf.Dir <- "/mnt/research/NMDL/MFM_WB/VariantCalling"
vcf.files <- list.files(vcf.Dir)[grep("vcf", list.files(vcf.Dir))]
genes <- unlist(lapply(strsplit(vcf.files, "_"), function(x) x[[1]][1]))
names(vcf.files) <- genes

#' Read VCF files into R
vcf <- lapply(genes, function(x) read.vcfR(paste(vcf.Dir, vcf.files[x], sep="/"), verbose = FALSE))
names(vcf) <- genes
vcf

#' Add ID column to VCF data by concating the chromosome and position of each cSNP
vcf <- lapply(vcf, function(x) addID(x, sep = "_"))
head(vcf[[1]])



#' ### Review Genotypes

#' Extract genotype matrix per gene
geno <- lapply(vcf, function(x) extract.gt(x, return.alleles=TRUE))

# Add animal names to column
anim <- unlist(lapply(strsplit(colnames(geno[[1]]), "[.]"), function(x) x[[1]][1]))
for (i in names(geno)){
    colnames(geno[[i]]) <- anim
}

#' Number of coding SNP per gene
unlist(lapply(geno, nrow))

#' ### Frequence of genotypes per gene
freq <- lapply(geno, function(x) apply(x, 1, table))


#' ### Save Genotype Information

#' Save called variant information to R data file
save(vcf.files, vcf, geno, freq, file=paste(getwd(), "gene_variants.Rdata", sep="/"))

#' Save frequency per gene

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

#' Write genotype matrix per gene
for (i in names(geno)){
    write.table(geno[[i]], file=paste(i, "genotypes.txt", sep="_"), 
        col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
}



#' ### Run R Script
#+ eval = FALSE
htmlRunR
Genotypes.R nodes=1,cpus-per-task=1,time=03:00:00,mem=2G \
+WB MFM Genotypes per Gene

