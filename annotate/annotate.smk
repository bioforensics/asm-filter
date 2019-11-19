import os

rule annotate:
    input:



rule prokka:
    input:
        '{sample}/analysis/asm_final.fasta'
    output:
        directory("{sample}/analysis/prokka/")
    shell:
        'prokka --outdir {output} --threads {threads} {input} '


rule abricate:
    input:
        '{sample}/analysis/asm_final.fasta'
    output:
        '{sample}/analysis/abricate/asm_final.abricate.out'
    shell:
        'abricate --db {config[abricate_db]} --threads {config[threads]} {input} > {output} '


        
