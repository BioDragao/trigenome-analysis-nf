#!/usr/bin/env nextflow

// /*
// ################
// NEXTFLOW Global Config
// ################
// */

// params.outdir = "results"

// /*
// ################
// GROOVY CODE (NATIVE)
// ################
// */

// name = Channel.from( 'Clojure', 'ClojureScript', 'Scheme', 'OCaml')

// process groovyPrintNames1 {
//     scratch true

//     input:
//     val name

//     exec:
//     println "Hello Mr. $name"
// }



// /*
// ################
// BASH SHELL
// ################
// */

// process bashShell {

//     output:
//     stdout bashShell_result

//     shell:
//     """
//     printf $SHELL
//     """
// }

// bashShell_result.println { it.trim() }

/*
################
gzip these files
################
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
gzip -dc G04868_R1.fastq.gz > G04868_R1.fastq

gzip -dc G04868_R2.fastq.gz > G04868_R2.fastq

    """
}

bashShell_result.println { it.trim() }




/*
###############
trimmomatic
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 G04868_1.fastq G04868_2.fastq G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

    """
}

bashShell_result.println { it.trim() }




/*
###############
bwa_index_reference_genome
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

bwa index NC000962_3.fasta

    """
}

bashShell_result.println { it.trim() }



/*
###############
map_and_generate_sam_file
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
bwa mem -R "@RG\tID:G04868\tSM:G04868\tPL:Illumina" -M NC000962_3.fasta G04868_1_trimmed_paired.fastq G04868_2_trimmed_paired.fastq > G04868.sam
    """
}

bashShell_result.println { it.trim() }



/*
###############
samtools_faidx_reference_genome
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
samtools faidx NC000962_3.fasta
    """
}

bashShell_result.println { it.trim() }







/*
###############
convert_sam_file_to_bam_file
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
samtools view -bt NC000962_3.fasta.fai G04868.sam > G04868.bam
    """
}

bashShell_result.println { it.trim() }



/*
###############
sort_bam_file
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
samtools sort G04868.bam -o G04868.sorted.bam
    """
}

bashShell_result.println { it.trim() }





/*
###############
samtools_index_sorted_bam
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
samtools index G04868.sorted.bam

    """
}

bashShell_result.println { it.trim() }





/*
###############
mapping_statistics
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

samtools flagstat G04868.sorted.bam > G04868_stats.txt

    """
}

bashShell_result.println { it.trim() }


/*
###############
samtools_mpileup
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
samtools mpileup -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta G04868.sorted.bam | bcftools call -O v -vm -o G04868.raw.vcf

    """
}

bashShell_result.println { it.trim() }




/*
###############
vcfutils_filter
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

vcfutils.pl varFilter -d 10 -D 2000 G04868.raw.vcf > G04868.filt.vcf

    """
}

bashShell_result.println { it.trim() }



/*
###############
bgzip_filt_file
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

bgzip -c G04868.filt.vcf > G04868.filt.vcf.gz

    """
}

bashShell_result.println { it.trim() }




/*
###############
run_tabix
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
tabix -p vcf G04868.filt.vcf.gz

    """
}

bashShell_result.println { it.trim() }




/*
###############
snpEff
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

java -Xmx4g -jar /opt/snpEff/snpEff.jar -no-downstream -no-upstream -v -c /opt/snpEff/snpEff.config NC000962_3 G04868.filt.vcf > G04868.ann.vcf.gz
    """
}

bashShell_result.println { it.trim() }




/*
###############
velveth_assembly
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velveth G04868_41 41 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq
    """
}

bashShell_result.println { it.trim() }




/*
###############
velvetg_produce_graph
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velvetg G04868_41 -exp_cov auto -cov_cutoff auto
    """
}

bashShell_result.println { it.trim() }





/*
###############
assemblathon_stats
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """
assemblathon_stats.pl ./G04868_41/contigs.fa
    """
}

bashShell_result.println { it.trim() }




/*
###############
velveth_assembly
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velveth G04868_49 49 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq

    """
}

bashShell_result.println { it.trim() }



/*
###############
velvetg_produce_graph
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velvetg G04868_49 -exp_cov auto -cov_cutoff auto

    """
}

bashShell_result.println { it.trim() }



/*
###############
assemblathon_stats
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

assemblathon_stats.pl ./G04868_49/contigs.fa

    """
}

bashShell_result.println { it.trim() }



/*
###############
velveth_assembly
###############
*/


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velveth G04868_55 55 -fastq -shortPaired  G04868_1_trimmed_paired.fastq G04868_1_trimmed_unpaired.fastq -fastq -short G04868_2_trimmed_paired.fastq G04868_2_trimmed_unpaired.fastq

    """
}

bashShell_result.println { it.trim() }


/*
###############
velvetg_produce_graph
###############
*/



process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

velvetg G04868_55 -exp_cov auto -cov_cutoff auto

    """
}

bashShell_result.println { it.trim() }



/*
###############
assemblathon_stats
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

assemblathon_stats.pl ./G04868_55/contigs.fa
    """

}

bashShell_result.println { it.trim() }


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

assemblathon_stats.pl ./G04868_41/contigs.fa
    """

}

bashShell_result.println {it.trim()}


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

assemblathon_stats.pl ./G04868_49/contigs.fa
    """

}

bashShell_result.println {it.trim()}


process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

assemblathon_stats.pl ./G04868_55/contigs.fa
    """

}

bashShell_result.println {it.trim()}

/*
###############
NOTE: Highest quality k_mer : 49
###############
*/


/*
###############
abacas_align_contigs
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

cd G04868_49 &&  cp ../NC000962_3.fasta ./ && abacas.pl -r ../NC000962_3.fasta -q contigs.fa -p promer -b -d -a
    """

}

bashShell_result.println {it.trim()}

/*
###############
prokka_annotation
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

cd ./G04868_49 && prokka --outdir ./G04868_prokka --prefix G04868 contigs.fa_NC000962_3.fasta.fasta
    """

}

bashShell_result.println {it.trim()}


/*
###############
gzip_compression
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

gzip -c G04868_1.fastq > G04868_1.fastq.gz
    """

}

bashShell_result.println {it.trim()}


/*
###############
gzip_compression
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

gzip -c G04868_2.fastq > G04868_2.fastq.gz
    """

}

bashShell_result.println {it.trim()}


/*
###############
snippy_command
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

snippy --cpus 4 --outdir G04868 --ref ./NC000962_3.gbk --R1 ./G04868_1.fastq.gz --R2 ./G04868_2.fastq.gz
    """

}

bashShell_result.println {it.trim()}


/*
###############
SNPtable
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

SNPtable_filter_Mtb.R core.tab
    """

}

bashShell_result.println {it.trim()}


/*
###############
HammingFasta
###############
*/

process bashShell {

    output:
    stdout bashShell_result

    shell:

    """

HammingFasta.R coreSNP_alignment_filtered.fas
    """

}

bashShell_result.println {it.trim()}

/*
###############
ALL DONE
###############
*/
