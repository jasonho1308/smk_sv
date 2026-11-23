# Opt-in DUX4 branches; independent of the normal multi-caller SV consensus.
from snakemake.exceptions import WorkflowError

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


rule dux4_dragen:
    input:
        expand("dux4_dragen/{sample}/attempt.json", sample=SAMPLES),


def dragen_dux4_bam(wildcards):
    pattern = config["dragen_dux4"].get("bam_pattern", rules.minimap2.output.bam)
    return pattern.replace("{sample}", wildcards.sample)


rule run_dragen_dux4:
    input:
        hash_table=config["dragen_dux4"]["hash_table"],
        bam=dragen_dux4_bam,
        runner=workflow.source_path("../scripts/run_dragen_dux4.py"),
    output:
        status="dux4_dragen/{sample}/attempt.json",
        results=directory("dux4_dragen/{sample}/dragen"),
    params:
        bypass="--skip-sanity-check" if config["dragen_dux4"].get("skip_sanity_check", True) else "",
    resources:
        dragen_slots=1,
    log:
        "logs/{sample}/dux4_dragen.log",
    shell:
        """
        python3 {input.runner:q} --sample {wildcards.sample:q} --bam {input.bam:q} \
            --hash-table {input.hash_table:q} --status {output.status:q} \
            --results {output.results:q} --log {log:q} {params.bypass}
        """
