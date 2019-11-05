
filePairs  = Channel.fromFilePairs('./*R{1,2}.fastq.gz')


//process echoFileNames {
//
//    echo true
//
//    input:
//    val file_list from filePairs
//
//    """
//   echo ${file_list[0]}
//   echo ${file_list[1][0]}
//   echo ${file_list[1][1]}
//   """
//}

process gzipFiles {

    echo true

    input:
    val fileList from filePairs


    script:
    oldName = fileList[1][0]
    newName = oldName.toString().split("\\.")[0] + ".fastq"

    """
    echo ${oldName}
    echo ${newName}
    """
}
