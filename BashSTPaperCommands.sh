#!/bin/bash

########################################################################
# Author: Anthony Hawkins                                              #
# Description: A collection of command line arguments, used to produce #
# various parts of the analysis.                                       #
########################################################################

#########################################
# Running Lace            ###############
#########################################

#Use assembled transcripts from Trinity (Trinity.fasta) and Clustering from Corset(D_inf-clusters-0.3.txt)
# Using 20 cores, and producing the annotated transcript annotation (--alternate) 
# Optionally giving the output directory
 
python Lace.py Trinity.fasta D_inf-clusters-0.3.txt --cores 20 --alternate --outputDir ST_20120709_trinity_corset


#########################################
# Making dynamic blocks   ###############
#########################################

#Concatenate the splice junctions for all samples (output of STAR) into one
#splice junction list
 
cat *.SJ.out.tab > SJ.out.tab

#Create dynamic blocks annotation based on splice junctions (SJ.out.tab) and superTranscript sequence (SuperDuper.fasta) 
python Mobius.py SJ.out.tab SuperDuper.fasta 


#########################################
# Running featureCounts   ###############
#########################################

# 20 threads (-T)
# Require both reads to aligne (-B)
# Allow reads which overlap multiple exons (-O)
# Use annotation file SuperDuper.gff (accepts .gff by defeault, but could also be .gtf or .saf), could be genome, or dynamic blocks instead
# Count at feature level (exon/blocks) rather than gene (-f)
# Output to Counts/counts.txt

featureCounts -T 20 -p -B -O --fraction  -a SuperDuper.gff  -f -o Counts/counts.txt *.out.bam


#########################################
# Running BLAT to find cluster gene #####
# correspondence for DTU analysis   #####
#########################################

blat SuperDuper.fasta HumanSuperDuper.fasta --minIdentity=98 output.psl


#########################################
# Running Salmon   ######################
#########################################

#Create salmon transcript index: transcripts.idx
salmon index -i transcripts.idx -t ../20120709_trinity/Trinity.fasta --type quasi -k 31

#Create a bash script to quantify transcript abundances per sample using Salmon and the previously created index
ls /mnt/storage/shared/public_data/cuffdiff2_data/SR*_1.trimmed.fastq  | rev | cut -c 17- | rev | uniq > myfiles
for i in `cat myfiles`; do
        echo "Performing Salmon quant on file $i"
        outname=${i##*/}
        salmon quant -i transcripts.idx -p 12 -l A -1 "$i"_1.trimmed.fastq -2 "$i"_2.trimmed.fastq -o "$outname"
done

#########################################
# Running Kallisto ######################
#########################################

#Create the kallisto transcript index: transcripts.idx
kallisto index -i transcripts.idx ../20120709_trinity/Trinity.fasta

#Create a bash script to quantify transcript abundances per sample using Kallisto and the previously created index
ls /mnt/storage/shared/public_data/cuffdiff2_data/SR*_1.trimmed.fastq  | rev | cut -c 17- | rev | uniq > myfiles
for i in `cat myfiles`; do
        echo "Performing kallisto quant on file $i"
        outname=${i##*/}
        kallisto quant -i transcripts.idx -t 8 -o "$outname" -b 100 "$i"_1.trimmed.fastq "$i"_2.trimmed.fastq
done

##########################################
# Running STAR  ##########################
##########################################


#Make Index for SuperTranscript
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir STARIndex --genomeFastaFiles SuperDuper.fasta --limitGenomeGenerateRAM 75070378709

#Get the files
firstpass=ST_20120709_trinity_corset/1pass
secondpass=ST_20120709_trinity_corset/2pass
GenomeDir=ST_20120709_trinity_corset/STARIndex
outdir=ST_20120709_trinity_corset/
GenomeFasta=ST_20120709_trinity_corset/SuperDuper.fasta
GenomeAnno=ST_20120709_trinity_corset/SuperDuper.gtf

#First pass mapping
cd /group/bioi1/shared/public_data/cuffdiff2_data/

STAR --genomeDir $GenomeDir \
--readFilesIn SRR493366_1.trimmed.fastq,SRR493367_1.trimmed.fastq,SRR493368_1.trimmed.fastq,SRR493369_1.trimmed.fastq,SRR493370_1.trimmed.fastq,SRR493371_1.trimmed.fastq  SRR493366_2.trimmed.fastq,SRR493367_2.trimmed.fastq,SRR493368_2.trimmed.fastq,SRR493369_2.trimmed.fastq,SRR493370_2.trimmed.fastq,SRR493371_2.trimmed.fastq \
--outFileNamePrefix $firstpass/SRR \
--outSAMtype BAM Unsorted --runThreadN 24;


#Generate genome with junctions from 1st pass
STAR --genomeDir $secondpass --runMode genomeGenerate --genomeFastaFiles $GenomeFasta --sjdbFileChrStartEnd $firstpass/SRRSJ.out.tab --sjdbOverhang 100 --runThreadN 20 --limitGenomeGenerateRAM 75070378709;

#Run 2nd pass mapping to the new genome defined from sjdb
ls *.fastq | rev | cut -c 17- | rev | uniq > myfiles
for i in `cat myfiles`; do
        echo "Performing alignment on file $i"
        STAR --genomeDir $secondpass \
        --readFilesIn "$i"_1.trimmed.fastq "$i"_2.trimmed.fastq \
        --outFileNamePrefix $outdir/"$i" \
        --outSAMtype BAM Unsorted --runThreadN 20 --outSJfilterOverhangMin 12 12 12 12;
done