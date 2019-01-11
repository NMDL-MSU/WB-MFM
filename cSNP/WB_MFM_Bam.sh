#==============================================================================
#   File: WB_MFM_Bam.sh
#   Directory code: /mnt/research/NMDL/MFM_WB
#   Date: 01/10/2019
#   Description: Copy BAM files from Warmblood MFM project to HPC and download 
#        reference cDNA sequence (FASTA file) from ensemble
#   Run: Interactively
#------------------------------------------------------------------------------
#   Input files in directory:
#       /media/velezdeb/D4E6-57FD/RNA-SEQ/analysis_data/WB_MFM_RNA_seq_2018/Analysis/STAR_run_uisng_equab2/*Aligned.sortedByCoord.out.bam
#
#   Output files to directory:
#       /mnt/research/NMDL/MFM_WB
#==============================================================================

### Copy BAM files to work directory from 8TB external hard drive
# NOTE: This code was run on the desktop computer terminal (not in HPC)

# Animal Index
cd /media/velezdeb/D4E6-57FD/RNA-SEQ/analysis_data/WB_MFM_RNA_seq_2018/Analysis/STAR_run_uisng_equab2
idx=(`ls *Aligned.sortedByCoord.out.bam | cut -f1 -dA`)

# Copy BAM to HPC
for ((i=0; i<${#idx[@]} ; i++ )) do
sshpass -f "/home/velezdeb/.pass" scp -p ${idx[$i]}Aligned.sortedByCoord.out.bam velezdeb@gateway.hpcc.msu.edu:/mnt/research/NMDL/MFM_WB
done

# Get reference used in generating BAM files (needed when call coding SNP)
sshpass -f "/home/velezdeb/.pass" scp -p horse_r86_all_corrected.fa* velezdeb@gateway.hpcc.msu.edu:/mnt/research/NMDL/MFM_WB

### Download EquCab2 FASTA for cDNA
cd /mnt/research/NMDL/MFM_WB
wget Index of ftp://ftp.ensembl.org/pub/release-94/fasta/equus_caballus/cdna/Equus_caballus.EquCab2.cdna.all.fa.gz
gunzip Equus_caballus.EquCab2.cdna.all.fa.gz

### Download chromosome lengths for EquCab2 (Assembly: GCA_000002305.1)
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/assembly/GCA_000/GCA_000002/GCA_000002305.1_sequence_report.txt

# Keep only information of chromosomed (remove unannotated scaffolds)
cat GCA_000002305.1_sequence_report.txt | cut -f1,2,3 | head -33 > EquCab2_chrom_length.txt

### Obtain genomic location of genes of interest
# Gene List
echo 'DES CRYAB MYOT PDLIM3 FLNC BAG3 ENSECAT00000027018.1 ENSECAG00000009005.1 PLEC LMNA ACTA1 HSPB8 KY PYROXD1 SQSTM1 TIA1' > genes.txt
sed -i 's/ /\n/g' genes.txt

# Extract genomic posirion from cDNA fasta file
idx=(`cat genes.txt`)
for ((i=0; i<${#idx[@]} ; i++ )) do
    grep -w ${idx[$i]} Equus_caballus.EquCab2.cdna.all.fa | cut -f1,3,7 -d' ' >> genes_location.txt
done

# Remove column names in each row and add to first row
sed -i 's/chromosome:EquCab2://' genes_location.txt 
sed -i 's/gene_symbol://' genes_location.txt 
sed -i '1i Ensembl Location Symbol' genes_location.txt 

# Add symbol to genes without gene symbol
sed -i 's/X:108319433:108323089:1/X:108319433:108323089:1 FHL1/' genes_location.txt
sed -i 's/4:107640843:107662925:1/4:107640843:107662925:1 DNAJB6/' genes_location.txt

