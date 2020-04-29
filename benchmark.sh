#!/bin/bash
for i in {001..050}
do
mkdir H3K27me3/E$i
mkdir H3K27ac/E$i
mkdir H3K9me3/E$i
#https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/E$i-H3K27me3.tagAlign.gz -O H3K27me3/E$i/input.tagAlign.gz
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/E$i-H3K27ac.tagAlign.gz -O H3K27ac/E$i/input.tagAlign.gz
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/E$i-H3K9me3.tagAlign.gz -O H3K9me3/E$i/input.tagAlign.gz

#https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E$i'_15_coreMarks_mnemonics.bed.gz' -O H3K27me3/E$i/15_coreMarks_mnemonics.bed.gz
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E$i'_15_coreMarks_mnemonics.bed.gz' -O H3K27ac/E$i/15_coreMarks_mnemonics.bed.gz
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E$i'_15_coreMarks_mnemonics.bed.gz' -O H3K9me3/E$i/15_coreMarks_mnemonics.bed.gz
done
for i in H3K27me3 H3K27ac H3K9me3; do for j in {001..050}; do mkdir $i/E$j/temp ; gzip -d $i/E$j/input.tagAlign.gz ; gzip -d $i/E$j/15_coreMarks_mnemonics.bed.gz ; done; done
for i in H3K27me3 H3K27ac H3K9me3; do for j in {001..050}; do snakemake $i/E$j/temp/OK --cores 8 ; done; done
