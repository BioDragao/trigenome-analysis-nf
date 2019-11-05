
filePairs  = Channel.fromFilePairs('./*R{1,2}.fastq.gz')


process gzipFiles {

    echo true

    input:
    val fileList from filePairs


    script:
    oldName = fileList[0]

    """
    echo ${oldName}
    """
}
