#!/bin/bash

for i in H3K27me3 H3K27ac H3K9me3; do for j in {001..050}; do mkdir $i/E$j/temp ; gzip -d $i/E$j/input.tagAlign.gz ; gzip -d $i/E$j/15_coreMarks_mnemonics.bed.gz ; done; done
for i in H3K27me3 H3K27ac H3K9me3; do for j in {001..050}; do snakemake $i/E$j/temp/OK --cores 8 ; done; done
