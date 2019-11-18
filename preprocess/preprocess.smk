import os
import snakemake
import fileinput

rule preprocess:
    input:
        expand('{sample}/analysis/mash_screen/results.tab', sample=config["samples"]),
        expand('{sample}/seq/filtered_{sample}_R{pair}.fastq.gz.cov{cov}.fastq', pair=config["pairs"], cov=config["covs"], sample=config["samples"])


rule mash_screen:
    input:
        read1='{sample}/seq/filtered_{sample}_R1.fastq.gz',
        read2='{sample}/seq/filtered_{sample}_R2.fastq.gz'
    output:
        '{sample}/analysis/mash_screen/results.tab'
    run:
        shell('mash screen -p {config[threads]} -w  {config[mashdb]} {input.read1} | sort -gr > {output}')

# create read sketch
# use mash output estimate genome size
# random downsample reads to coverages given in config file
rule size_downsample:
    input:
        expand('{{sample}}/seq/filtered_{{sample}}_R{pair}.fastq.gz', pair=config["pairs"])
    output:
        expand('{{sample}}/seq/filtered_{{sample}}_R{pair}.fastq.gz.cov{{cov}}.fastq', pair=config["pairs"]),
        sketch='{sample}/analysis/mash_sketch/read_sketch.{cov}.msh'
    run:
        if config[genome_size] == 0:
            for line in shell("genome_estimation.py {input} {output.sketch}", iterable=True):
                genome_size = line
        else:
            genome_size = config[genome_size]
        shell("sample-reads {input} {genome_size} {wildcards.cov}")
        

# legacy support for cutadapt & fastqc
# prefer fastp be run instead see below
# run cutadapt if set                
if "cutadapt" == config["qa"]:
    rule trim:
        input:
            read1='{sample}/seq/{sample}_R1.fastq.gz',
            read2='{sample}/seq/{sample}_R2.fastq.gz'
        output:
            read1='{sample}/seq/filtered_{sample}_R1.fastq.gz',
            read2='{sample}/seq/filtered_{sample}_R2.fastq.gz',
	    report='{sample}/seq/qa.html'
        log:         
            '{sample}/seq/qa.log'
        run:
            shell('cutadapt -a XXX -A XXX -m {minlength} -q {config[qual]},{config[qual]} -o {output.read1} -p {output.read1} {input.read1} {input.read2} > {log} 2>&1 ')
            shell('fastqc -o {sample}/seq {input.read1} {input.read2}')

else:
    # run fastp by default
    # default qual 20
    # rewrite html file to find javascript on GDNA
    rule trim:
        input:
            read1='{sample}/seq/{sample}_R1.fastq.gz',
            read2='{sample}/seq/{sample}_R2.fastq.gz'
        output:
            read1='{sample}/seq/filtered_{sample}_R1.fastq.gz',
            read2='{sample}/seq/filtered_{sample}_R2.fastq.gz',
	    report='{sample}/seq/qa.html'
        log:         
            '{sample}/seq/qa.log'
        # fastp max threads is 8
        threads: 8
        run:
            shell('fastp -w {threads} -l {config[minlength]} -q {config[qual]} --in1 {input.read1} --in2 {input.read2} --out1 {output.read1} --out2 {output.read2} --html {output.report} -j {wildcards.sample}/seq/fastp.json > {log} 2>&1')
            shell("sed -i 's|https://cdn.plot.ly/plotly-latest.min.js|file:///nbacc/data/plotly-latest.min.js|g' {output.report}")



            


        
