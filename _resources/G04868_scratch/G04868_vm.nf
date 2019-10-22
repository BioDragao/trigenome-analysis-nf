
filePairs  = Channel.fromFilePairs('./*R{1,2}.fastq.gz')

process bar {

    echo true

    input:
    val file_list from filePairs

   """
   echo ${file_list[0]}
   echo ${file_list[1][0]}
   echo ${file_list[1][1]}
   """
}

