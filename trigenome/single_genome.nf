//====== genome files needed ============
// NOTE

/*
G04868_L003_R1.fastq.gz
G04868_L003_R2.fastq.gz
NC000962_3.fasta
NC000962_3.gbk
*/

// TODO rely on the input param to have global variabel of genomeName

//====== gzip decompress ============
// DONE

/*
gzippedFilePairsCh = Channel.fromFilePairs('G04868_L003_R{1,2}.fastq.gz')


process gzipDecompressFiles {

    echo true

    input:

    val fileList from gzippedFilePairsCh

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
*/


//====== trimmomatic ============
// DONE

/*
fastqFilePairsCh = Channel.fromFilePairs('G04868_L003_R{1,2}.fastq.gz')


process trimmomatic {

    echo true

    input:

    val fileList from fastqFilePairsCh


    script:

    for (int i = 0; i < fileList.size(); i++) {
        oldR1Name = fileList[i + 1][0]
        newR1Paired = oldR1Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
        newR1Unpaired = oldR1Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"

        oldR2Name = fileList[i + 1][1]
        newR2Paired = oldR2Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
        newR2Unpaired = oldR2Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"


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
*/

//====== bwa index ============
// DONE
/*

referenceGenomeCh = Channel.fromPath("./NC000962_3.fasta")

process bwaIndexReferenceGenome {
//    conda 'bwa'
//    conda './tese.yaml'

    echo true

    input:
    val refGenome from referenceGenomeCh

    script:

    """
    bwa index ${refGenome}
    """
}

*/
////======== samtools_faidx_reference_genome =======
//// DONE
//
//
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
//    samtools faidx  ${refGenome}
//    """
//}
//
//

//======= map_and_generate_sam_file =======
// DONE

fastqFilePairsCh = Channel.fromFilePairs('G04868_L003_R{1,2}_trimmed_paired.fastq')
referenceGenomeCh = Channel.fromPath("./NC000962_3.fasta")


process mapAndGenerateSamFile {
//    conda 'bwa'
//    conda './tese.yaml'

    echo true

    input:
    val refGenome from referenceGenomeCh
    val fastqFiles from fastqFilePairsCh

    script:

    samFileName = fastqFiles[1][0].toString().split("\\.")[0].split("\\_")[0] + "_" + fastqFiles[1][0].toString().split("\\.")[0].split("\\_")[1] + ".sam"
    fastqPairedFile1 = fastqFiles[1][0]
    fastqPairedFile2 = fastqFiles[1][1]

    """
    bwa mem -R "@RG\\tID:G04868\\tSM:G04868\\tPL:Illumina" -M ${refGenome} ${fastqPairedFile1} ${fastqPairedFile2} > ${
        samFileName
    }
    """
}


////======== convert_sam_file_to_bam_file =======
//
//genomeFromPathCh = Channel.fromPath('./G04868_L003_R1.fastq.gz')
//referenceGenomeFaiCh = Channel.fromPath('./NC000962_3.fasta.fai')
//samFileCh = Channel.fromPath("./G04868_L003.sam")
//
//process convertSamFileToBamFile {
////    conda 'bwa'
////    conda './tese.yaml'
//
//    echo true
//
//    input:
//    val genomeFromPath from genomeFromPathCh
//    val referenceGenomeFai from referenceGenomeFaiCh
//    val samFile from samFileCh
//
//
//    script:
//
//    genomeName = "G04868_" + genomeFromPath.toString().split("\\.")[0].split("\\_")[1]
//    bamFile = "./" + genomeName + ".bam"
//
//    """
//    samtools view -bt ${referenceGenomeFai} ${samFile} > ${bamFile}
//    """
//}
//
////======== sort_bam_file =======
////// TODO
////
////
////
////
////
////
////======== samtools_index_sorted_bam =======
////// TODO
////
////
////
////
////
////
////======== mapping_statistics =======
////// TODO
////
////
////
////
////
////
////======== samtools_mpileup =======
////// TODO
////
////
////
////
////
////
////======== vcfutils_filter =======
////// TODO
////
////
////
////
////
////======== bgzip_filt_file =======
////// TODO
////
////
////
////
////
////
////======== run_tabix =======
////// TODO
////
////
////
////
////
////
////======== snpEff =======
////// TODO
////
////
////
////
////======== velveth_assembly =======
////// TODO
////
////
////
////
////
////======== velvetg_produce_graph =======
////// TODO
////
////
////
////
////
////
////======== assemblathon_stats =======
////// TODO
////
////
////
////
////
////
////======== velveth_assembly =======
////// TODO
////
////
////
////
////
////
////======== velvetg_produce_graph =======
////// TODO
////
////
////
////
////
////======== assemblathon_stats =======
////// TODO
////
////
////
////
////
////======== velveth_assembly =======
////// TODO
////
////
////
////======== velvetg_produce_graph =======
////// TODO
////
////
////
////
////
////======== assemblathon_stats =======
////// TODO
////
////
////
////
////======== abacas_align_contigs =======
////// TODO
////
////
////
////
////======== prokka_annotation =======
////// TODO
////
////
////
////
////
////======== gzip_compression =======
////// TODO
////
////
////
////
////
////======== snippy_command =======
////// TODO
////
////
////
////#========================================================
////# The next section starts when we've done the analysis for all genomes
////#========================================================
////
////======== snippy_core =======
////
////
////======== SNPtable =======
////// TODO
////
////
////
////
////======== HammingFasta =======
////// TODO
////
////
