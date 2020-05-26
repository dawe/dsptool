#Run by:
# $ snakemake -s /pathtosnakefile/Snakefile justdoit --cores 12 --config path=/pathtoinputfolderwithsubfolders

import glob

def find_tagAlign(path_of_tagAligns):
    mylist=glob.glob('{}/**/*.tagAlign'.format(path_of_tagAligns), recursive=True)
    print ("mylist is:", mylist)
    NAME, PDIR = [], []
    for element in mylist:
        name=element.rstrip(".tagAlign").split("/")[-1]
        parent=element.rstrip(".tagAlign").rstrip('/'+name)
        NAME.append(name)
        PDIR.append(parent)
    return NAME,PDIR

NAME, PDIR = find_tagAlign(config['path'])
# print ("NAME IS:",NAME)
# print ("PDIR IS:",PDIR)

rule r1:
    input: "{pdir}/{name}.tagAlign"
    output: "{pdir}/{name}.bam"
    shell: "bedToBam -i {input[0]} -g hg19.chrom.sizes | samtools sort -o {output[0]}"

rule r2:
    input: "{pdir}/{name}.bam"
    output: "{pdir}/{name}.bam.bai"
    shell: "samtools index {input[0]}"

rule r3:
    input: "{pdir}/{name}.bam", "{pdir}/{name}.bam.bai"
    output: "{pdir}/{name}.bw"
    shell: "bamCoverage -b {input[0]} -o {output[0]} --extend 300 "

rule r4:
    input: "{pdir}/{name}.bw"
    output: "{pdir}/{name}_denoised.bw", "{pdir}/temp/{name}_denoised_highresolution.bed", "{pdir}/temp/{name}_denoised_lowresolution.bed"
    shell: "python dsp.py -i {input[0]} -d {output[0]} -temp {wildcards.pdir}/temp"

rule r5:
    input: "{pdir}/temp/{name}_denoised_highresolution.bed" , "{pdir}/temp/{name}_denoised_lowresolution.bed"
    output: "{pdir}/{name}_segmented.bed", "{pdir}/temp/{name}_temp.bed"
    shell:
        """
        bedtools intersect -a {input[1]} -b {input[0]} -wa -wb | bedtools groupby -c 4,5,6 -o distinct,min,max | cut -f 4- > {output[1]}
        bedtools intersect -a {output[1]} -b {input[0]} -wa -wb | bedtools groupby -c 5,6 -o collapse,collapse > {output[0]}
        """

rule justdoit:
    input: expand("{pdir}/{name}_segmented.bed", zip, pdir=PDIR, name=NAME)
