filePairs = Channel.fromFilePairs('./*_R{1,2}.fastq.gz')


process gzipFiles {

    echo true

    input:

    val fileList from filePairs


    output:

    file "${fileList[1][0].toString().split(" \\.")[0]}.fastq" into unzippedFiles


    script:

    """
    gzip -dc ${fileList[1][0]} > ${fileList[1][0].toString().split("\\.")[0]}.fastq

    gzip -dc ${fileList[1][1]} > ${fileList[1][1].toString().split("\\.")[0]}.fastq

    """
}


unzippedFiles
        .flatMap()
        .subscribe { println "${it.name}" }


//
//process splitLetters {
//
//    echo true
//
//    input:
//
//    val file from unzippedFiles
//
//    script:
//
//    """
//    echo ${file}
//    """
//}


