#!/usr/bin/env bash

rm single_genome_report.html
rm single_genome_timeline.html
rm single_genome_flowchart_dag.html
rm .nextflow.*

nextflow single_genome.nf -with-report single_genome_report.html -with-timeline single_genome_timeline.html -with-dag single_genome_flowchart_dag.html # -with-trace

git add -f ./single_genome.sh ./single_genome.nf \
#./single_genome_report.html  \
#./single_genome_timeline.html  \
#./single_genome_flowchart_dag.html \
./.nextflow.log
