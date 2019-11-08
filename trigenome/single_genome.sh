#!/usr/bin/env bash

rm single_genome.html
rm .nextflow.*


nextflow  single_genome.nf  -with-report single_genome.html

git add -f ./single_genome.sh ./single_genome.nf ./single_genome.html ./.nextflow.log
