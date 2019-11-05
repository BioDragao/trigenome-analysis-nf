//====== globals ============

referenceGenome = Channel.fromPath("./NC000962_3.fasta")


//====== trimmomatic ============

/*
gzippedFilePairs = Channel.fromFilePairs('./*_R{1,2}.fastq.gz')


process gzipFiles {

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

*/

//====== trimmomatic ============

//fastqFilePairs = Channel.fromFilePairs('./*_R{1,2}.fastq')
//
//
//process trimmomatic {
//
//    echo true
//
//    input:
//
//    val fileList from fastqFilePairs
//
//    // TODO implement the <output> to pass the unzipped files to a channel in the pipeline
//
//    script:
//
//    for (int i = 0; i < fileList.size(); i++) {
//        oldR1Name = fileList[i + 1][0]
//        newR1Paired = oldR1Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
//        newR1Unpaired = oldR1Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"
//
//        oldR2Name = fileList[i + 1][1]
//        newR2Paired = oldR2Name.toString().split("\\.")[0] + "_trimmed_paired.fastq"
//        newR2Unpaired = oldR2Name.toString().split("\\.")[0] + "_trimmed_unpaired.fastq"
//
//
////  java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 G04868_1.fastq G04868_2.fastq G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
//
//        return """
//            echo java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE  -phred33  \
//            ${oldR1Name}  \
//            ${oldR2Name}  \
//            ${newR1Paired}  \
//            ${newR1Unpaired}  \
//            ${newR2Paired}  \
//            ${newR2Unpaired} \
//            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
//
//            """
//    }
//}

//====== bwa index ============


process bwaIndexReferenceGenome {
//    conda 'bwa'
//    conda '/some/path/my-env.yaml'

    echo true

    input:
    val refGenome from referenceGenome

    script:

    """
    echo bwa index ${refGenome}
    """
}

