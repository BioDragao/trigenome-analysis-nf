#!/usr/bin/env bash
git add -f *nf *sh *log *html


rm single_genome.html
rm .nextflow.*
nextflow  single_genome.nf  -with-report single_genome.html

