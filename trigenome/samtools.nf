//======== samtools_faidx_reference_genome =======

// TODO this is repeated, get value without consuming the content of channel in above usage
referenceGenome = Channel.fromPath("./NC000962_3.fasta")

process samtoolsFaidxReferenceGenome {
//    conda 'bwa'
//    conda './tese.yaml'

    echo true

    input:
    val refGenome from referenceGenome

    script:

//    samtools faidx NC000962_3.fasta

    """
    echo samtools faidx  ${refGenome}
    """
}

//
//
//======== convert_sam_file_to_bam_file =======

referenceGenome = Channel.fromPath("./NC000962_3.fasta")

process samtoolsFaidxReferenceGenome {
//    conda 'bwa'
//    conda './tese.yaml'

    echo true

    input:
    val refGenome from referenceGenome

    script:

//        samtools view -bt NC000962_3.fasta.fai G04868.sam > G04868.bam

    """
    echo samtools faidx  ${refGenome}
    """
}

//
//
//
//======== sort_bam_file =======
//
//
//        samtools sort G04868.bam -o G04868.sorted.bam
//
//
//
//
//======== samtools_index_sorted_bam =======
//
//
//        samtools index G04868.sorted.bam
//
//
//
//
//======== mapping_statistics =======
//
//
//        samtools flagstat G04868.sorted.bam > G04868_stats.txt
//
//
//
//
//======== samtools_mpileup =======
//
//
//        samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta G04868.sorted.bam | bcftools call -O v -vm -o G04868.raw.vcf
//
//
//
//