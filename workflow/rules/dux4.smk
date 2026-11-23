# Opt-in screen; independent of the normal multi-caller SV consensus.

rule dux4:
    input:
        expand("dux4/{sample}/{sample}.junctions.tsv", sample=SAMPLES),


rule scan_dux4:
    conda:
        "../envs/dux4.yaml"
    input:
        bam=rules.minimap2.output.bam,
        index=rules.minimap2.output.bai,
        reference=config["fasta"],
        fai=rules.faidx.output.fai,
        partners=workflow.source_path("../resources/dux4/chm13v2.0.partners.bed"),
        sources=workflow.source_path("../resources/dux4/chm13v2.0.sources.bed"),
        scan=workflow.source_path("../scripts/scan_dux4.py"),
    output:
        junctions="dux4/{sample}/{sample}.junctions.tsv",
        reads="dux4/{sample}/{sample}.reads.tsv",
        loci="dux4/{sample}/{sample}.loci.tsv",
    threads: 4
    resources:
        mem_mb=2048,
    log:
        "logs/{sample}/dux4.log",
    shell:
        """
        python {input.scan:q} --require-chm13 \
            --bam {input.bam:q} --bam-index {input.index:q} --reference {input.reference:q} \
            --partners {input.partners:q} --sources {input.sources:q} \
            --junctions {output.junctions:q} --reads {output.reads:q} --loci {output.loci:q} \
            --threads {threads} > {log:q} 2>&1
        """
