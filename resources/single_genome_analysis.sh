#!/usr/bin/env bash

# gzip decompress <<<<<
# DONE
gzip -dc G04868_L003_R1.fastq.gz >G04868_L003_R1.fastq

gzip -dc G04868_L003_R2.fastq.gz >G04868_L003_R2.fastq

# trimmomatic <<<<<
# DONE

java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 G04868_L003_R1.fastq G04868_L003_R2.fastq G04868_L003_R1_trimmed_paired.fastq G04868_L003_R1_trimmed_unpaired.fastq G04868_L003_R2_trimmed_paired.fastq G04868_L003_R2_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# bwa_index_reference_genome <<<<<
# DONE

bwa index NC000962_3.fasta

# samtools_faidx_reference_genome <<<<<
# DONE

samtools faidx NC000962_3.fasta

# map_and_generate_sam_file <<<<<
# DONE

bwa mem -R "@RG\tID:G04868\tSM:G04868\tPL:Illumina" -M NC000962_3.fasta G04868_L003_R1_trimmed_paired.fastq G04868_L003_R2_trimmed_paired.fastq > G04868_L003.sam

# convert_sam_file_to_bam_file <<<<<

samtools view -bt NC000962_3.fasta.fai G04868_L003.sam > G04868_L003.bam

# sort_bam_file <<<<<

samtools sort G04868_L003.bam -o G04868_L003.sorted.bam

# samtools_index_sorted_bam <<<<<

samtools index G04868_L003.sorted.bam

# mapping_statistics <<<<<

samtools flagstat G04868_L003.sorted.bam >G04868_L003_stats.txt

# samtools_mpileup <<<<<

samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta G04868_L003.sorted.bam | bcftools call -O v -vm -o G04868_L003.raw.vcf

# vcfutils_filter<<<<<

vcfutils.pl varFilter -d 10 -D 2000 G04868_L003.raw.vcf >G04868_L003.filt.vcf

# bgzip_filt_file <<<<<

bgzip -c G04868_L003.filt.vcf >G04868_L003.filt.vcf.gz

# run_tabix <<<<<

tabix -p vcf G04868_L003.filt.vcf.gz

# snpEff <<<<<

java -Xmx4g -jar /opt/snpEff/snpEff.jar -no-downstream -no-upstream -v -c /opt/snpEff/snpEff.config NC000962_3 G04868_L003.filt.vcf >G04868_L003.ann.vcf.gz

# velveth_assembly <<<<<

velveth G04868_L003_41 41 -fastq -shortPaired G04868_L003_R1_trimmed_paired.fastq G04868_L003_R1_trimmed_unpaired.fastq -fastq -short G04868_L003_R2_trimmed_paired.fastq G04868_L003_R2_trimmed_unpaired.fastq

# velvetg_produce_graph <<<<<

velvetg G04868_L003_41 -exp_cov auto -cov_cutoff auto

# assemblathon_stats <<<<<

assemblathon_stats.pl ./G04868_L003_41/contigs.fa

# velveth_assembly <<<<<

velveth G04868_L003_49 49 -fastq -shortPaired G04868_L003_R1_trimmed_paired.fastq G04868_L003_R1_trimmed_unpaired.fastq -fastq -short G04868_L003_R2_trimmed_paired.fastq G04868_L003_R2_trimmed_unpaired.fastq

# velvetg_produce_graph <<<<<

velvetg G04868_L003_49 -exp_cov auto -cov_cutoff auto

# assemblathon_stats <<<<<

assemblathon_stats.pl ./G04868_L003_49/contigs.fa

# velveth_assembly <<<<<

velveth G04868_L003_55 55 -fastq -shortPaired G04868_L003_R1_trimmed_paired.fastq G04868_L003_R1_trimmed_unpaired.fastq -fastq -short G04868_L003_R2_trimmed_paired.fastq G04868_L003_R2_trimmed_unpaired.fastq

# velvetg_produce_graph <<<<<

velvetg G04868_L003_55 -exp_cov auto -cov_cutoff auto

# assemblathon_stats <<<<<

assemblathon_stats.pl ./G04868_L003_55/contigs.fa

assemblathon_stats.pl ./G04868_L003_41/contigs.fa

assemblathon_stats.pl ./G04868_L003_49/contigs.fa

assemblathon_stats.pl ./G04868_L003_55/contigs.fa

# Highest quality k_mer : 49

# abacas_align_contigs <<<<<

cd G04868_L003_49 && cp ../NC000962_3.fasta ./ && abacas.pl -r ../NC000962_3.fasta -q contigs.fa -p promer -b -d -a

# prokka_annotation <<<<<

cd ./G04868_L003_49 && prokka --outdir ./G04868_L003_prokka --prefix G04868_L003 contigs.fa_NC000962_3.fasta.fasta

# gzip_compression <<<<<

gzip -c G04868_L003_R1.fastq > G04868_L003_R1.fastq.gz

# gzip_compression <<<<<

gzip -c G04868_L003_R2.fastq > G04868_L003_R2.fastq.gz

# snippy_command <<<<<

snippy --cpus 4 --outdir G04868_L003 --ref ./NC000962_3.gbk --R1 ./G04868_L003_R1.fastq.gz --R2 ./G04868_L003_R2.fastq.gz

#========================================================
# The next section starts when we've done the analysis for all genomes
#========================================================

# snippy_core <<<<<

snippy-core G04868_L003

# SNPtable <<<<<

SNPtable_filter_Mtb.R core.tab

# HammingFasta <<<<<

HammingFasta.R coreSNP_alignment_filtered.fas

# ALL DONE
