#!/usr/bin/env bash

rm .nextflow.*
nextflow  single_genome.nf  -with-report single_genome.html

