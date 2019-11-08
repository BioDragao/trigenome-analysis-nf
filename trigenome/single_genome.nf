//====== globals ============



//====== gzip ============
// DONE

gzippedFilePairs = Channel.fromFilePairs('G04868_L003_R{1,2}.fastq.gz')


process gzipDecompressFiles {

    echo true

    input:

    val fileList from gzippedFilePairs

    // TODO implement the <output> to pass the unzipped files to a channel in the pipeline

    script:

    for (int i = 0; i < fileList.size(); i++) {
        oldR1Name = fileList[i + 1][0]
        newR1Name = oldR1Name.toString().split("\\.")[0]

        oldR2Name = fileList[i + 1][1]
        newR2Name = oldR2Name.toString().split("\\.")[0]

        return """
            gzip -dc ${oldR1Name} > ${newR1Name}.fastq

            gzip -dc ${oldR2Name} > ${newR2Name}.fastq
            """
    }
}


//====== trimmomatic ============
// DONE
fastqFilePairs = Channel.fromFilePairs('G04868_L003_R{1,2}.fastq.gz')


process trimmomatic {

    echo true

    input:

    val fileList from fastqFilePairs

    // TODO implement the <output> to pass the unzipped files to a channel in the pipeline

    script:

    for (int i = 0; i < fileList.size(); i++) {
        oldR1Name = fileList[i + 1][0]
        newR1Paired = oldR1Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
        newR1Unpaired = oldR1Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"

        oldR2Name = fileList[i + 1][1]
        newR2Paired = oldR2Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
        newR2Unpaired = oldR2Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"


//  java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 G04868_1.fastq G04868_2.fastq G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

        return """
            java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE  -phred33  \
            ${oldR1Name}  \
            ${oldR2Name}  \
            ${newR1Paired}  \
            ${newR1Unpaired}  \
            ${newR2Paired}  \
            ${newR2Unpaired} \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

            """
    }
}

//====== bwa index ============
// DONE

referenceGenome = Channel.fromPath("./NC000962_3.fasta")

process bwaIndexReferenceGenome {
//    conda 'bwa'
//    conda './tese.yaml'

    echo true

    input:
    val refGenome from referenceGenome

    script:

    """
    bwa index ${refGenome}
    """
}

//======= map_and_generate_sam_file =======
//
//        bwa mem - R "@RG\tID:G04868\tSM:G04868\tPL:Illumina" -M NC000962_3 . fasta G04868_1_trimmed_paired.fastq G04868_2_trimmed_paired . fastq > G04868.sam
//
//
//======== samtools_faidx_reference_genome =======
//
//
//        samtools faidx NC000962_3 . fasta

// TODO this is repeated, get value without consuming the content of channel in above usage
//referenceGenome = Channel.fromPath("./NC000962_3.fasta")
//
//process samtoolsFaidxReferenceGenome {
////    conda 'bwa'
////    conda './tese.yaml'
//
//    echo true
//
//    input:
//    val refGenome from referenceGenome
//
//    script:
//
//    """
//    echo samtools faidx  ${refGenome}
//    """
//}
//
//
//
//
//======== convert_sam_file_to_bam_file =======
//
//
//        samtools view - bt NC000962_3 . fasta.fai G04868.sam > G04868.bam
//
//
//
//
//======== sort_bam_file =======
//
//
//        samtools sort G04868 . bam - o G04868.sorted.bam
//
//
//
//
//======== samtools_index_sorted_bam =======
//
//
//        samtools index G04868 . sorted.bam
//
//
//
//
//======== mapping_statistics =======
//
//
//        samtools flagstat G04868 . sorted.bam > G04868_stats.txt
//
//
//
//
//======== samtools_mpileup =======
//
//
//        samtools mpileup - Q 23 -d 2000 -C 50 -ugf NC000962_3 . fasta G04868.sorted.bam | bcftools call -O v -vm - o G04868 . raw.vcf
//
//
//
//
//======== vcfutils_filter =======
//
//
//        vcfutils.pl varFilter - d 10 -D 2000 G04868.raw.vcf > G04868.filt.vcf
//
//
//
//
//======== bgzip_filt_file =======
//
//
//        bgzip -c G04868.filt.vcf > G04868.filt.vcf.gz
//
//
//
//
//======== run_tabix =======
//
//
//        tabix -p vcf G04868 . filt.vcf.gz
//
//
//
//
//======== snpEff =======
//
//
//        java -Xmx4g -jar / opt / snpEff / snpEff.jar -no- downstream - no - upstream - v - c / opt / snpEff / snpEff.config NC000962_3 G04868 . filt.vcf > G04868.ann.vcf.gz
//
//
//======== velveth_assembly =======
//
//
//        velveth G04868_41 41 -fastq -shortPaired G04868_1_trimmed_paired . fastq G04868_1_trimmed_unpaired.fastq - fastq - short G04868_2_trimmed_paired . fastq G04868_2_trimmed_unpaired.fastq
//
//
//
//======== velvetg_produce_graph =======
//
//
//        velvetg G04868_41 -exp_cov auto -cov_cutoff auto
//
//
//
//
//======== assemblathon_stats =======
//
//
//        assemblathon_stats.pl./G04868_41/ contigs.fa
//
//
//
//
//======== velveth_assembly =======
//
//
//        velveth G04868_49 49 -fastq -shortPaired G04868_1_trimmed_paired . fastq G04868_1_trimmed_unpaired.fastq - fastq - short G04868_2_trimmed_paired . fastq G04868_2_trimmed_unpaired.fastq
//
//
//
//
//======== velvetg_produce_graph =======
//
//
//        velvetg G04868_49 -exp_cov auto -cov_cutoff auto
//
//
//
//
//======== assemblathon_stats =======
//
//
//        assemblathon_stats.pl./G04868_49/ contigs.fa
//
//
//
//======== velveth_assembly =======
//
//
//        velveth G04868_55 55 -fastq - shortPaired G04868_1_trimmed_paired . fastq G04868_1_trimmed_unpaired.fastq - fastq - short G04868_2_trimmed_paired . fastq G04868_2_trimmed_unpaired.fastq
//
//======== velvetg_produce_graph =======
//
//
//        velvetg G04868_55 - exp_cov auto -cov_cutoff auto
//
//
//
//======== assemblathon_stats =======
//
//
//        assemblathon_stats.pl./G04868_55/ contigs.fa
//
//
//assemblathon_stats.pl./G04868_41/ contigs.fa
//
//
//assemblathon_stats.pl./G04868_49/ contigs.fa
//
//
//assemblathon_stats.pl./G04868_55/ contigs.fa
//
//
//#
//Highest quality
//k_mer: 49
//
//
//======== abacas_align_contigs =======
//
//
//        cd G04868_49 && cp../NC000962_3.fasta ./ && abacas.pl - r../NC000962_3.fasta -q contigs.fa -p promer -b -d -a
//
//
//======== prokka_annotation =======
//
//
//cd ./ G04868_49 && prokka-- outdir./G04868_prokka --prefix G04868 contigs.fa_NC000962_3.fasta.fasta
//
//
//
//======== gzip_compression =======
//
//
//gzip -c G04868_1.fastq > G04868_1.fastq.gz
//
//
//
//======== gzip_compression =======
//
//
//gzip -c G04868_2.fastq > G04868_2.fastq.gz
//
//
//
//======== snippy_command =======
//
//
//snippy --cpus 4 --outdir G04868 --ref ./ NC000962_3 . gbk-- R1./G04868_1.fastq.gz --R2 ./ G04868_2 . fastq.gz
//
//
//
//======== SNPtable =======
//
//
//        SNPtable_filter_Mtb.R core.tab
//
//
//
//======== HammingFasta =======
//
//
//        HammingFasta.R coreSNP_alignment_filtered.fas
