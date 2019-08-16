#!/usr/bin/env nextflow


/*
################
NOTE

Before running this, ensure
 - virtualbox or bioconda setup
 - rclone setup
 - nextflow and JDK-8 (possibly graalVM)
################
*/




/*
################
NEXTFLOW Global Config
################
*/

params.outdir = "results"

/*
################
GROOVY CODE (NATIVE)
################
*/

name = Channel.from('Emilyn', 'Abhinav' )

process groovyCode {
    scratch true

    input:
    val name

    exec:
    println "$name"
}



/*
################
gzip these files
################
*/

process runGzip {

    output:
    stdout runGzip_result

    shell:

    """
gzip -dc G04868_R1.fastq.gz > G04868_R1.fastq

gzip -dc G04868_R2.fastq.gz > G04868_R2.fastq

    """
}

runGzip_result.println { it.trim() }



/*
###############
trimmomatic
###############
*/


process runTrimmomatic {

    output:
    stdout runTrimmomatic_result

    shell:

    """

java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 G04868_1.fastq G04868_2.fastq G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

    """
}

runTrimmomatic_result.println { it.trim() }


// /*
// ###############
// bwa_index_reference_genome
// ###############
// */

// process runBwaIndex {

//     output:
//     stdout runBwaIndex_result

//     shell:

//     """

// bwa index NC000962_3.fasta

//     """
// }

// runBwaIndex_result.println { it.trim() }



// /*
// ###############
// map_and_generate_sam_file
// ###############
// */



// process runBwaMapping {

//     output:
//     stdout runBwaMapping_result

//     shell:

//     """
// bwa mem -R "@RG\tID:G04868\tSM:G04868\tPL:Illumina" -M NC000962_3.fasta G04868_1_trimmed_paired.fastq G04868_2_trimmed_paired.fastq > G04868.sam
//     """
// }

// runBwaMapping_result.println { it.trim() }


// /*
// ###############
// samtools_faidx_reference_genome
// ###############
// */


// process runSamtoolsFaidx {

//     output:
//     stdout runSamtoolsFaidx_result

//     shell:

//     """
// samtools faidx NC000962_3.fasta
//     """
// }

// runSamtoolsFaidx_result.println { it.trim() }


// /*
// ###############
// convert_sam_file_to_bam_file
// ###############
// */



// process runSamToBam {

//     output:
//     stdout runSamToBam_result

//     shell:

//     """
// samtools view -bt NC000962_3.fasta.fai G04868.sam > G04868.bam
//     """
// }

// runSamToBam_result.println { it.trim() }



// /*
// ###############
// sort_bam_file
// ###############
// */


// process runSortBam {

//     output:
//     stdout runSortBam_result

//     shell:

//     """
// samtools sort G04868.bam -o G04868.sorted.bam
//     """
// }

// runSortBam_result.println { it.trim() }


// /*
// ###############
// samtools_index_sorted_bam
// ###############
// */


// process runSamtoolsIndex {

//     output:
//     stdout runSamtoolsIndex_result

//     shell:

//     """
// samtools index G04868.sorted.bam

//     """
// }

// runSamtoolsIndex_result.println { it.trim() }





// /*
// ###############
// mapping_statistics
// ###############
// */

// process runSamtoolsStat {

//     output:
//     stdout runSamtoolsStat_result

//     shell:

//     """

// samtools flagstat G04868.sorted.bam > G04868_stats.txt

//     """
// }

// runSamtoolsStat_result.println { it.trim() }


// /*
// ###############
// samtools_mpileup
// ###############
// */


// process runSamtoolsMpileup {

//     output:
//     stdout runSamtoolsMpileup_result

//     shell:

//     """
// samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta G04868.sorted.bam | bcftools call -O v -vm -o G04868.raw.vcf

//     """
// }

// runSamtoolsMpileup_result.println { it.trim() }

// /*
// ###############
// vcfutils_filter
// ###############
// */


// process runVcfutils {

//     output:
//     stdout runVcfutils_result

//     shell:

//     """

// vcfutils.pl varFilter -d 10 -D 2000 G04868.raw.vcf > G04868.filt.vcf

//     """
// }

// runVcfutils_result.println { it.trim() }



// /*
// ###############
// bgzip_filt_file
// ###############
// */



// process runBgzip {

//     output:
//     stdout runBgzip_result

//     shell:

//     """

// bgzip -c G04868.filt.vcf > G04868.filt.vcf.gz

//     """
// }

// runBgzip_result.println { it.trim() }

// /*
// ###############
// run_tabix
// ###############
// */

// process runTabix {

//     output:
//     stdout runTabix_result

//     shell:

//     """
// tabix -p vcf G04868.filt.vcf.gz

//     """
// }

// runTabix_result.println { it.trim() }


// /*
// ###############
// snpEff
// ###############
// */



// process runSnpeff {

//     output:
//     stdout runSnpeff_result

//     shell:

//     """

// java -Xmx4g -jar /opt/snpEff/snpEff.jar -no-downstream -no-upstream -v -c /opt/snpEff/snpEff.config NC000962_3 G04868.filt.vcf > G04868.ann.vcf.gz
//     """
// }

// runSnpeff_result.println { it.trim() }

// /*
// ###############
// velveth_assembly
// ###############
// */

// process runVelvetH41 {

//     output:
//     stdout runVelvetH41_result

//     shell:

//     """

// velveth G04868_41 41 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq
//     """
// }

// runVelvetH41_result.println { it.trim() }


// /*
// ###############
// velvetg_produce_graph
// ###############
// */


// process runVelvetG41 {

//     output:
//     stdout runVelvetG41_result

//     shell:

//     """

// velvetg G04868_41 -exp_cov auto -cov_cutoff auto
//     """
// }

// runVelvetG41_result.println { it.trim() }

// /*
// ###############
// assemblathon_stats
// ###############
// */


// process runAssemblathon41 {

//     output:
//     stdout runAssemblathon41_result

//     shell:

//     """
// assemblathon_stats.pl ./G04868_41/contigs.fa
//     """
// }

// runAssemblathon41_result.println { it.trim() }


// /*
// ###############
// velveth_assembly
// ###############
// */

// process runVelvetH49 {

//     output:
//     stdout runVelvetH49_result

//     shell:

//     """

// velveth G04868_49 49 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq

//     """
// }

// runVelvetH49_result.println { it.trim() }


// /*
// ###############
// velvetg_produce_graph
// ###############
// */



// process runVelvetG49 {

//     output:
//     stdout bashShell_result

//     shell:

//     """

// velvetg G04868_49 -exp_cov auto -cov_cutoff auto

//     """
// }

// runVelvetG49_result.println { it.trim() }



// /*
// ###############
// assemblathon_stats
// ###############
// */


// process runAssemblathon49 {

//     output:
//     stdout runAssemblathon49_result

//     shell:

//     """

// assemblathon_stats.pl ./G04868_49/contigs.fa

//     """
// }

// runAssemblathon49_result.println { it.trim() }


// /*
// ###############
// velveth_assembly
// ###############
// */


// process runVelvetH55 {

//     output:
//     stdout runVelvetH55_result

//     shell:

//     """

// velveth G04868_55 55 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq

//     """
// }

// runVelvetH55_result.println { it.trim() }


// /*
// ###############
// velvetg_produce_graph
// ###############
// */



// process runVelvetG55 {

//     output:
//     stdout runVelvetG55_result

//     shell:

//     """

// velvetg G04868_55 -exp_cov auto -cov_cutoff auto

//     """
// }

// runVelvetG55_result.println { it.trim() }




// /*
// ###############
// assemblathon_stats
// ###############
// */


// process runAssemblathon55 {

//     output:
//     stdout runAssemblathon55_result

//     shell:

//     """

// assemblathon_stats.pl ./G04868_50/contigs.fa

//     """
// }

// runAssemblathon55_result.println { it.trim() }

// /*
// ###############
// Find the highest quality k-mer using the assemblathon-results
// ###############
// */


// /*
// NOTE: Highest quality k_mer : 49
// */

// assemblathon_stats.pl ./G04868_41/contigs.fa
// assemblathon_stats.pl ./G04868_49/contigs.fa
// assemblathon_stats.pl ./G04868_55/contigs.fa


// /*
// ###############
// abacas_align_contigs
// ###############
// */

// process runAbacas {

//     output:
//     stdout runAbacas_result

//     shell:

//     """

// cd G04868_49 &&  cp ../NC000962_3.fasta ./ && abacas.pl -r ../NC000962_3.fasta -q contigs.fa -p promer -b -d -a
//     """

// }

// runAbacas_result.println {it.trim()}

// /*
// ###############
// prokka_annotation
// ###############
// */

// process runProkka {

//     output:
//     stdout runProkka_result

//     shell:

//     """

// cd ./G04868_49 && prokka --outdir ./G04868_prokka --prefix G04868 contigs.fa_NC000962_3.fasta.fasta
//     """

// }

// runProkka_result.println {it.trim()}


// /*
// ###############
// gzip_compression
// ###############
// */

// process runGzipCompression {

//     output:
//     stdout runGzipCompression_result

//     shell:

//     """

// gzip -c G04868_1.fastq > G04868_1.fastq.gz

// gzip -c G04868_2.fastq > G04868_2.fastq.gz
//     """

// }

// runGzipCompression_result.println {it.trim()}



// /*
// ###############
// snippy_command
// ###############
// */

// process runSnippy {

//     output:
//     stdout runSnippy_result

//     shell:

//     """
// snippy --cpus 4 --outdir G04868 --ref ./NC000962_3.gbk --R1 ./G04868_1.fastq.gz --R2 ./G04868_2.fastq.gz
//     """

// }

// runSnippy_result.println {it.trim()}


// /*
// ###############
// SNPtable
// ###############
// */

// process runSnpTable {

//     output:
//     stdout runSnpTable_result

//     shell:

//     """

// SNPtable_filter_Mtb.R core.tab
//     """

// }

// runSnpTable_result.println {it.trim()}


// /*
// ###############
// HammingFasta
// ###############
// */

// process runHammingFasta {

//     output:
//     stdout runHammingFasta_result

//     shell:

//     """

// HammingFasta.R coreSNP_alignment_filtered.fas
//     """

// }

// runHammingFasta_result.println {it.trim()}

// /*
// ###############
// ALL DONE
// ###############
// */
