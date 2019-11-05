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
