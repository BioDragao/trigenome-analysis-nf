#!/usr/bin/env nextflow

//================================
// General Algorithms
// - read fromFilePairs in the channel
// -
//
//
//
//
//
//================================





// /*
// ################
// NEXTFLOW Global Config
// ################
// */

 params.outdir = "G04868_analysis"
 params.genomeName = "G04868"


// =========================

process sayHello {

    "echo 'Hello, World!"
}


// /*
// ################
// Read file pairs
// ################
// */



// gzippedFiles = Channel
//     .fromFilePairs('../_resources/G04868_scratch/*_R{1,2}.fastq.gz')

// process echoFileNames {

//     input:
//     tuple val(x) from gzippedFiles

//     "echo $x"
// }


/*
################
echo file names
################
*/


// gzipped = Channel.fromPath('../_resources/G04868_scratch/*.fastq.gz')

// process echoFileNames {

//     input:
//     val x from gzipped

//     exec:
//     println "$x"
// }




process runGzip {

    output:
    stdout runGzip_result

    shell:

    """
gzip -dc G04868_R1.fastq.gz > G04868_R1.fastq

gzip -dc G04868_R2.fastq.gz > G04868_R2.fastq

    """
}



