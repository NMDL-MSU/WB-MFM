#==============================================================================
#   File: call_gene_variants.sh
#   Directory code: /mnt/research/NMDL/MFM_WB/VariantCalling
#   Date: 01/10/2019
#   Description: Run mpileup and bcftools on bam files on select genes for 
#          Warmblood MFM project 
#   Run: bash call_gene_variants.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       /mnt/research/NMDL/MFM_WB
#
#       Files:
#           *Aligned.sortedByCoord.out.bam
#           horse_r86_all_corrected.fa
#           genes_location.txt
#
#   Output files to directory:
#       /mnt/research/NMDL/MFM_WB/VariantCalling
#==============================================================================

# Create work directory
cd /mnt/research/NMDL/MFM_WB
mkdir VariantCalling
dir=/mnt/research/NMDL/MFM_WB/VariantCalling
mv call_gene_variants.sh $dir

# Create qstat directory
mkdir VariantCalling/qstat
qstat=/mnt/research/NMDL/MFM_WB/VariantCalling/qstat

# Bam files directory
Bam=/mnt/research/NMDL/MFM_WB

# Rename Bam files (use only animal ID which transfer to VCF id names) and index with samtools
module purge
module load powertools
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SAMtools/1.9

idx=(`ls *.bam`)
id=(`ls *.bam | cut -f1 -d_`)
for ((i=0; i<${#idx[@]} ; i++ )) do
    mv ${idx[$i]} ${id[$i]}.bam 
    samtools index ${id[$i]}.bam
done


# Reference EquCab2 (same as the reference used for generating BAM files)
ref=/mnt/research/NMDL/MFM_WB/horse_r86_all_corrected.fa

# Gene indexes
tail -23 genes_location.txt | cut -f2 -d' ' | cut -f1,2 -d: > start
tail -23 genes_location.txt | cut -f2 -d' ' | cut -f3 -d: > end
paste start end | sed 's/\t/-/' > index_genes.txt
rm start end
index=(`cat index_genes.txt`)
nms=(`tail -23 genes_location.txt | cut -f3 -d' '`)

# Submit script to HPC
for ((i=0; i<${#nms[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=50G
#SBATCH -J '${nms[$i]}'
#SBATCH -o '${nms[$i]}.o%j'

# Load Modules
module purge
module load powertools
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SAMtools/1.9
module load bcftools/1.9.64
#module load tabix/0.2.6
module load matplotlib/2.1.2-Python-3.6.4
module list

# Move to output directory
cd '$Bam'

# Run mpileup
bcftools mpileup -Ou -C50 -E -Q25 -a DV,AD,ADF,ADR,SP -f '$ref' -r '${index[$i]}' '*.bam' \
    | bcftools call -mv -Oz \
    | bcftools filter -s LowQual -g3 -G10 \
    -e '"'"'%QUAL<30'"'"' \
    -o '$dir'/'${nms[$i]}'_'${index[$i]}'.vcf

# Job details
echo 'Job Details'
scontrol show job $SLURM_JOB_ID' > $qstat/${nms[$i]}_${index[$i]}.qsub

# Submit script to hpcc
cd $qstat
sbatch ${nms[$i]}_${index[$i]}.qsub

done

