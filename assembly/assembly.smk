import os
import shutil

rule assembly:
    input:
        expand('{sample}/analysis/asm_final.fasta', sample=config["samples"]),
        expand('{sample}/analysis/circle/circles.fasta', sample=config["samples"]),
        expand('{sample}/analysis/nn/nn.fasta', sample=config["samples"])

#calculate nearest nieghbor
rule nn:
    input:
        '{sample}/analysis/asm_final.fasta'
    output:
        '{sample}/analysis/nn/nn.fasta'
    run:
        out_dir = os.path.dirname(output[0])
        shell("mash dist -p {config[threads]} {config[mashdb_asm]} {input} | sort -gk 3 > {out_dir}/nn.list")
        shell("cp $(head -1 {out_dir}/nn.list | cut -f1 ) {output}")


rule circle:
    input:
        '{sample}/analysis/asm_final.fasta'
    output:
        '{sample}/analysis/circle/circles.fasta'
    run:
        out_dir = os.path.dirname(output[0])
        shell("circlep.py {input} {out_dir} > {output}")

if "pilon" in config["polish"]:        
   rule pilon:
       input:
           asm='{sample}/analysis/{assembler}.{cov}/asm.fasta',
           bam='{sample}/analysis/{assembler}.{cov}/asm.fasta.bam'
       output:
           '{sample}/analysis/asm_final.fasta'
       params:
           '{sample}/analysis/{assembler}.{cov}/pilon'
       shell:
           'pilon --genome {input.asm} --frags {input.bam} --outdir {params}'

rule scoring:
     input:
         expand('{{sample}}/analysis/{assembler}.{cov}/score', cov=config["covs"], assembler=config["assemblers"])
     output:
         '{sample}/analysis/asm_final.fasta'
     run:
         shell("score.sh {wildcards.sample}")
        
rule coverage_eval_filter:
     input:
         asm='{sample}/analysis/{assembler}.{cov}/asm.fasta',
         read1='{sample}/seq/filtered_{sample}_R1.fastq.gz.cov{cov}.fastq',
         read2='{sample}/seq/filtered_{sample}_R2.fastq.gz.cov{cov}.fastq'
     params:
         path="{sample}/analysis/{assembler}.{cov}/"
     output:
         '{sample}/analysis/{assembler}.{cov}/score',
         '{sample}/analysis/{assembler}.{cov}/asm.fasta.bam'
         
     shell:
         'coverage_eval.sh {input} {params.path} {config[threads]}'

         
if "spades" in config["assemblers"]:
   rule spades:
       input:
           read1='{sample}/seq/filtered_{sample}_R1.fastq.gz.cov{cov}.fastq',
           read2='{sample}/seq/filtered_{sample}_R2.fastq.gz.cov{cov}.fastq'
       output:
           '{sample}/analysis/spades.{cov}/asm.fasta'
       message: "Running spades "
       run:
           out_dir = os.path.dirname(output[0])

           shell("spades.py --careful -t {config[threads]} -m {config[memory]} -1 {input.read1} -2 {input.read2} {config[spades_options]} --cov-cutoff auto -o {out_dir}")
           shutil.copyfile(out_dir + "/scaffolds.fasta", output[0])

if "megahit" in config["assemblers"]:
    rule megahit:
        input:
            read1='{sample}/seq/filtered_{sample}_R1.fastq.gz.cov{cov}.fastq',
            read2='{sample}/seq/filtered_{sample}_R2.fastq.gz.cov{cov}.fastq'
        output:
            '{sample}/analysis/spades.{cov}/asm.fasta'
        message: "Running megahit "
        run:
            out_dir = os.path.dirname(output[0])

            shell("megahit -t {config[threads]} -1 {input.read1} -2 {input.read2} {config[megahit_options]} -o {out_dir}")
            shutil.copyfile(out_dir + "/final.contigs.fasta", output[0])


if "unicycler" in config["assemblers"]:
    rule unicycler:
        input:
            read1='{sample}/seq/filtered_{sample}_R1.fastq.gz.cov{cov}.fastq',
            read2='{sample}/seq/filtered_{sample}_R2.fastq.gz.cov{cov}.fastq'
        output:
            '{sample}/analysis/spades.{cov}/asm.fasta'
        message: "Running unicycler "        
        run:
            out_dir = os.path.dirname(output[0])

            shell("unicycler -t {config[threads]} -1 {input.read1} -2 {input.read2} {config[megahit_options]} -o {out_dir}")
            shutil.copyfile(out_dir + "/assembly.fasta", output[0])
        
