
filePairs  = Channel.fromFilePairs('./*R{1,2}.fastq.gz')


process gzipFiles {

    echo true

    input:
    val fileList from filePairs

    script:
    firstFileOldName = fileList[1][0]
    firstFileNewName = firstFileOldName.toString().split("\\.")[0] + ".fastq"

    secondFileOldName = fileList[1][1]
    secondFileNewName = secondFileOldName.toString().split("\\.")[0] + ".fastq"

    """
   gzip -dc ${firstFileOldName}   >  ${firstFileNewName}
   gzip -dc ${secondFileOldName}  >  ${secondFileNewName}
    """
}
