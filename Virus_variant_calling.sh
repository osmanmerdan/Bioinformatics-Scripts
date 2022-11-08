#!/bin/bash

#Author: Osman Merdan
#Date Last Modified: 08/11/2022

# Stop on any error 
# -e: Exit immediately if a command exits with a non-zero status.
# -u  Treat unset variables as an error when substituting.
# -x  Print commands and their arguments as they are executed.
set -ue




#########################
#   Preparing Files     #
#########################

# Enter the accession number for the reference genome
# -r option: do not allow backslashes to escape any characters
read -r -p "Input the refrence genome GenBank Accession Number:" AAC
read -r -p "Input SRR accession numbers as one column txt file:" SRR

# Make the directory for reference 
mkdir -p "$HOME"/"${AAC}"Project/refs

# The GenBank refrence file 
GBK="$HOME/${AAC}Project/refs/${AAC}.gb"
# The FASTA sequence 

REF="$HOME/${AAC}Project/refs/${AAC}.fa"

# Fetch the sequence 
bio fetch "$AAC"> "$GBK"

# Transform GenBank into FASTA 
bio fasta "$GBK"> "$REF"

# Build an index for the genome so that we can view in IGV
samtools faidx "$REF"

# Build the index
bwa index "$REF"

# Create a directory for raw reads and go into it  
#NC_045512
mkdir -p "$HOME"/"${AAC}"Project/raw

# Get the sequences which are listed in the srrlist.txt from SRR archive. 
echo "******** Raw Seqs Download Step Started **********"
parallel --bar "fastq-dump -X 20000 --split-files --outdir ~/${AAC}Project/raw {}" < "$SRR"
echo "******** Raw Seqs Downloaded **********"

# Fastqc 
mkdir ~/"${AAC}"Project/qc
parallel "fastqc ~/${AAC}Project/raw/{}_1.fastq ~/${AAC}Project/raw/{}_2.fastq -o ~/${AAC}Project/qc" < "$SRR"

# Trimming 
while read -r i;
do
trimmomatic PE\
 ~/"${AAC}"Project/raw/"${i}"_1.fastq\
 ~/"${AAC}"Project/raw/"${i}"_2.fastq\
 ~/"${AAC}"Project/raw/"${i}"_paired_1.fastq\
 ~/"${AAC}"Project/raw/"${i}"_unpaired_1.fastq\
 ~/"${AAC}"Project/raw/"${i}"_paired_2.fastq\
 ~/"${AAC}"Project/raw/"${i}"_unpaired_2.fastq\
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:7:25 MINLEN:150
 # :2 = seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
 # Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches.
 # :30 = palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment. 
 # :10 = simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
 # :2 = minAdapterLength: In addition to the alignment score, palindrome mode can verify that a minimum length of adapter has been detected. 
 # If unspecified, this defaults to 8 bases, for historical reasons. 
 # However, since palindrome mode has a very low false positive rate, this can be safely reduced, even down to 1, to allow shorter adapter fragments to be removed.
 # windowSize:7
 # requiredQuality:25
done < "$SRR"




#########################
#      Alignment        #
#########################

# Alignment with bwa mem (6 thread) 
while read -r i; do
echo "******* Mapping of the $i is started $(date +%d-%m-%Y_%H-%M-%S) *******"
bwa mem -t 6\
 "$REF"\
 ~/"${AAC}"Project/raw/"${i}"_paired_1.fastq\
 ~/"${AAC}"Project/raw/"${i}"_paired_2.fastq > ~/"${AAC}"Project/"${i}".sam
echo "****** Mapping of the $i is completed $(date +%d-%m-%Y_%H-%M-%S) ******"
done <"$SRR"

# Alignment filtering
# 1) samtool fixmate -r Remove unmapped reads and secondary alignments
# 2) sort alignments according to chr coordinates
# 3) filter min-map quality 30, exclude 2048 flag (supplementary alignment), select only proper pairs flag 2 
parallel "samtools fixmate -O bam -r ~/${AAC}Project/{}.sam - |\
 samtools sort - |\
 samtools view -q 30 -F 2048 -f 2 -b -h -O bam -o ~/${AAC}Project/{}.bam" :::: "$SRR"

# Alignment matrics 
# samtools stats -c: Coverage distribution min,max,step [1,1000,1]
parallel "samtools flagstat ~/${AAC}Project/{}.bam > ~/${AAC}Project/qc/{}-flagstat.txt" :::: "$SRR"
parallel "samtools stats -c 1,30000,50 ~/${AAC}Project/{}.bam > ~/${AAC}Project/qc/{}.stats" :::: "$SRR"
# Extracting coverage data: grep ^COV | cut -f 2- 

# Duplicate removal (if required)
#picard MarkDuplicates --INPUT SRR13710914_sorted.bam\
# --OUTPUT SRR13710914_deduplicated.bam\
# --REMOVE-DUPLICATES true




#########################
#    Variant Calling    #
#########################

# Indexing bam files 
parallel --bar "samtools index ~/${AAC}Project/{}.bam" :::: "$SRR"

# Variant calling and filtering 
# bcftools call --keep-alts: Keep all possible alternate alleles at variant sites
# bcftools call -m Alternative model for multiallelic and rare-variant calling (conflicts with -c)
# bcftools norm: Left-align and normalize indels; check if REF alleles match the reference;
# split multiallelic sites into multiple rows; recover multiallelics from multiple rows.
# bcftools norm -m : Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]
# bcftools norm -f : Refrence file 
# bcftools filter select variants with %90 of the reads supporting alternate allel

parallel --bar "bcftools mpileup --output-type v --max-depth 2000 --fasta-ref $REF ~/${AAC}Project/{}.bam |\
 bcftools call --ploidy 1 --variants-only --keep-alts -m --output-type v |\
 bcftools norm -m - -f $REF -|\
 bcftools filter -i '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]))>0.9 && %QUAL>30 && DP>100' > ~/${AAC}Project/{}.vcf" :::: "$SRR"



#########################
#   Variant Annotation  #
#########################

# Downloading snpEff database
# First line of fasta file indicates version number of the genome
version_number="$(grep -o '[A-Z]\+.[0-9]\+.[0-9]' "$REF")"
snpEff download "$version_number"

# Annotate variants
# To simplyfy output use '-no-downstream -no-upstream'
parallel "snpEff eff $version_number {} > {.}_annotated.vcf" ::: ~/"${AAC}"Project/*.vcf

# Extracting fields using SnpSift into tabular format for downstream analysis
# To extract fields one per line use the script from https://github.com/pcingola/SnpEff/blob/master/scripts/vcfEffOnePerLine.pl 
# Further documentation here https://pcingola.github.io/SnpEff/ss_extractfields/
wget https://raw.githubusercontent.com/pcingola/SnpEff/master/scripts/vcfEffOnePerLine.pl -P ~/"${AAC}"Project/

parallel --bar "cat {}\
|perl  ~/${AAC}Project/vcfEffOnePerLine.pl\
|SnpSift extractFields - -s ',' -e '.' 'CHROM'\
 'POS' 'ID' 'REF' 'ALT' 'FILTER' 'ANN[*].EFFECT'\
  'ANN[*].IMPACT' 'ANN[*].GENE' 'ANN[*].GENEID'\
   'ANN[*].FEATURE' 'ANN[*].FEATUREID'\
    'ANN[*].HGVS_C' 'ANN[*].HGVS_P' > {.}.txt" ::: ~/"${AAC}"Project/*_annotated.vcf


