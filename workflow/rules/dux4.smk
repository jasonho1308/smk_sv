# Opt-in DUX4 branches; independent of the normal multi-caller SV consensus.
include: "dux4_local_assembly.smk"


rule dux4:
    input:
        expand("dux4/{sample}/{sample}.junctions.tsv", sample=SAMPLES),
        expand("dux4/{sample}/{sample}.blast.tsv", sample=SAMPLES),
        expand("dux4/{sample}/{sample}.blast_filtered.tsv", sample=SAMPLES),


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
        blast="dux4/{sample}/{sample}.blast.tsv",
        filtered_blast="dux4/{sample}/{sample}.blast_filtered.tsv",
    threads: 16
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
            --blast {output.blast:q} --filtered-blast {output.filtered_blast:q} \
            --threads {threads} > {log:q} 2>&1
        """


# Manual target: run both Conda branches for one suspected sample.
rule dux4_followup:
    input:
        rules.scan_dux4.output,
        rules.local_assemble_dux4.output,
        rules.align_contigs_to_ref.output,
    output:
        manifest="dux4/{sample}/dux4_files.txt",
    shell:
        "printf '%s\\n' {input:q} > {output.manifest:q}"
