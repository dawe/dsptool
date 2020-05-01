import sys
if any([x for x in sys.argv[1].split("/") if x == 'H3K27me3']):
    features='11_BivFlnk\|10_TssBiv\|12_EnhBiv\|13_ReprPC\|14_ReprPCWk'
elif any([x for x in sys.argv[1].split("/") if x == 'H3K9me3']):
    features='9_Het\|15_Quies\|8_ZNF/Rpts'
else:
    features='1_TssA\|2_TssAFlnk\|6_EnhG\|7_Enh'


rule toBW:
    input:
        "{IN}.tagAlign"
    output:
        A="{IN}.bw",
        B="{IN}.bam"
    shell:
        """
        bedToBam -i {input} -g hg19.chrom.sizes | samtools sort -o {output.B}
        samtools index {output.B}
        bamCoverage -b {output.B} -o {output.A} --extend 300
        """

rule DSP:
    input:
        "{IN}input.bw"
    output:
        A="{IN}output.bw",
        B="{IN}temp/Wazowski.bed"
    shell:
        """
        python dsp.py -i {input} -d {output.A} -temp {wildcards.IN}/temp
        head {output.B}
        mkdir {wildcards.IN}temp/plt
        """

rule mnemonics:
    input:
        "{IN}/15_coreMarks_mnemonics.bed"
    output:
        "{IN}/temp/mnemonics.bed"
    shell:
        """
        # gzip -d {input}.gz
        grep '{features}' {input} > {output}
        grep '11_BivFlnk' {input} > {wildcards.IN}/temp/BivFlnk.bed
        grep '10_TssBiv' {input} > {wildcards.IN}/temp/TssBiv.bed
        grep '12_EnhBiv' {input} > {wildcards.IN}/temp/EnhBiv.bed
        grep '13_ReprPC' {input} > {wildcards.IN}/temp/ReprPC.bed
        grep '14_ReprPCWk' {input} > {wildcards.IN}/temp/ReprPCWk.bed
        """

rule GMM:
    input:
        A="{DIR}/temp/mnemonics.bed",
        B="{DIR}/temp/Wazowski.bed"
    output:
        C="{DIR}/temp/over-gmm-mnem.bed",
        D="{DIR}/temp/GMM.bed"
    shell:
        """
        python gmm.py {wildcards.DIR}/input.bw {input.B} {output.D}
        bedtools intersect -a {wildcards.DIR}/temp/Salivan.bed -b {output.D} -wa -wb | bedtools groupby -c 4,5,6 -o distinct,min,max | cut -f 4- > {wildcards.DIR}/temp/Lizard.bed
        bedtools intersect -a {wildcards.DIR}/temp/Lizard.bed -b {output.D} -wa -wb | bedtools groupby -c 5,6 -o collapse,collapse > {wildcards.DIR}/temp/GMMresult.bed
        bedtools intersect -wao -a {wildcards.DIR}/temp/GMMresult.bed -b {input.A} | bedtools groupby -o sum -c 10 > {output.C}
        rm {wildcards.DIR}/temp/Lizard.bed
        """


rule intDSP:
    input:
        A="{DIR}/mnemonics.bed",
        B="{DIR}/Wazowski.bed"
    output:
        C="{DIR}/over-dsp-mnem.bed",
        D="{DIR}/Segmented.bed"
    shell:
        """
        bedtools intersect -a {wildcards.DIR}/Salivan.bed -b {input.B} -wa -wb | bedtools groupby -c 4,5,6 -o distinct,min,max | cut -f 4- > {wildcards.DIR}/Lizard.bed
        bedtools intersect -a {wildcards.DIR}/Lizard.bed -b {input.B} -wa -wb | bedtools groupby -c 5,6 -o collapse,collapse > {output.D}
        bedtools intersect -wao -a {output.D} -b {input.A} | bedtools groupby -o sum -c 10 > {output.C}
        rm {wildcards.DIR}/Lizard.bed
        """

rule DShuffle:
    input:
        A="{DIR}/Segmented.bed",
        B="{DIR}/mnemonics.bed"
    output:
        "{DIR}/over-d-random-mnem.bed",
    shell:
        """
        bedtools shuffle -i {input.A} -g hg19.genome -noOverlapping | bedtools sort > {wildcards.DIR}/DShuffled.bed
        bedtools intersect -wao -a {wildcards.DIR}/DShuffled.bed -b {input.B} | bedtools groupby -o sum -c 10 > {output}
        """

rule SShuffle:
    input:
        A="{DIR}/sicer.bed",
        B="{DIR}/mnemonics.bed"
    output:
        "{DIR}/over-s-random-mnem.bed",
    shell:
        """
        bedtools shuffle -i {input.A} -g hg19.genome -noOverlapping | bedtools sort > {wildcards.DIR}/SShuffled.bed
        bedtools intersect -wao -a {wildcards.DIR}/SShuffled.bed -b {input.B} | bedtools groupby -o sum -c 8 > {output}
        """

rule GShuffle:
    input:
        A="{DIR}/GMM.bed",
        B="{DIR}/mnemonics.bed"
    output:
        "{DIR}/over-g-random-mnem.bed",
    shell:
        """
        bedtools shuffle -i {input.A} -g hg19.genome -noOverlapping | bedtools sort > {wildcards.DIR}/GShuffled.bed
        bedtools intersect -wao -a {wildcards.DIR}/GShuffled.bed -b {input.B} | bedtools groupby -o sum -c 8 > {output}
        """

rule MShuffle:
    input:
        "{DIR}/mnemonics.bed"
    output:
        "{DIR}/over-m-random-mnem.bed"
    shell:
        """
        bedtools shuffle -i {input} -g hg19.genome -noOverlapping | bedtools sort > {wildcards.DIR}/MShuffled.bed
        bedtools intersect -wao -a {wildcards.DIR}/MShuffled.bed -b {input} | bedtools groupby -o sum -c 9 > {output}
        """

rule Sicer:
    input:
        "{DIR}/input.bam"
    output:
        "{DIR}/temp/sicer.bed"
    shell:
        """
        sicer -t {input} -s hg19 -w 5000 -g 50000 -o {wildcards.DIR}
        cut -f 1-3 {wildcards.DIR}/input-W5000-G50000.scoreisland > {output}
        rm {wildcards.DIR}/input.bed
        """

rule intSicer:
    input:
        A="{DIR}/sicer.bed",
        B="{DIR}/mnemonics.bed"
    output:
        "{DIR}/over-sicer-mnem.bed"
    shell:
        """
        bedtools merge -i {input.A} > {wildcards.DIR}/x.bed
        bedtools intersect -wao -a {wildcards.DIR}/x.bed -b {input.B} | bedtools groupby -o sum -c 8 > {output}
        rm {wildcards.DIR}/x.bed
        """

# plot

rule overHMM:
    input:
        A="{DIR}/over-sicer-mnem.bed",
        B="{DIR}/over-m-random-mnem.bed",
        C="{DIR}/over-s-random-mnem.bed",
        D="{DIR}/over-g-random-mnem.bed",
        E="{DIR}/over-d-random-mnem.bed",
        F="{DIR}/over-dsp-mnem.bed",
        G="{DIR}/over-gmm-mnem.bed"
    output:
        Z="{DIR}/osm.txt",
        Y="{DIR}/omrm.txt",
        X="{DIR}/osrm.txt",
        W="{DIR}/ogrm.txt",
        V="{DIR}/odrm.txt",
        U="{DIR}/odm.txt",
        T="{DIR}/ogm.txt"
    shell:
        """
        awk "{{ print \$4/(\$3 - \$2) }}" {input.A} > {output.Z}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.B} > {output.Y}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.C} > {output.X}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.D} > {output.W}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.E} > {output.V}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.F} > {output.U}
        awk "{{ print \$4/(\$3 - \$2) }}" {input.G} > {output.T}
        """


rule plotOXM:
    input:
        S="{X}/osm.txt",
        R="{X}/omrm.txt",
        G="{X}/ogm.txt",
        D="{X}/odm.txt",
        RS="{X}/osrm.txt",
        RG="{X}/ogrm.txt",
        RD="{X}/odrm.txt"
    output:
        '{X}/plt/oxm.png'
    run:
        AXTITLE='\#,basepair,\(last,col\),\/,length,of,consider,region,\(col3-col2\),in,intersect,bed,files'
        TITLE="intersect,-wao,-a,Segmented.bed,-b,mnemonics.bed"
        COLNAME='Sicer,ShMnem,GMM,DSP,ShSicer,ShGMM,ShDSP'
        shell("""
        python boxplot.py {input.S} {input.R} {input.G} {input.D} {input.RS} {input.RG} {input.RD} {output} {TITLE} {AXTITLE} {COLNAME}
        """)

rule runall2:
    input:
        '{X}/plt/oxm.png'
    output:
        "{X}/OK"
    shell:
        """
        ls {wildcards.X} > {output}
        """
