# Manual DUX4 follow-up: request dux4/local_asm/<sample>/breakpoints.tsv.
DUX4_ASM = config["dux4_local_assembly"]
DUX4_ASM_DIR = "dux4/local_asm/{sample}"


def dux4_asm_bam(wildcards):
    pattern = DUX4_ASM.get("bam_pattern", rules.minimap2.output.bam)
    return pattern.replace("{sample}", wildcards.sample)


def dux4_asm_bam_index(wildcards):
    if "bam_pattern" not in DUX4_ASM:
        return rules.minimap2.output.bai.replace("{sample}", wildcards.sample)
    bam = dux4_asm_bam(wildcards)
    return DUX4_ASM.get("bam_index_pattern", bam + ".bai").replace("{sample}", wildcards.sample)


rule extract_dux4_igh_reads:
    # Historical rule name: collects reads at all configured partners, not just IGH.
    input:
        bam=dux4_asm_bam,
        index=dux4_asm_bam_index,
        fai=rules.faidx.output.fai,
        sources=workflow.source_path("../resources/dux4/chm13v2.0.assembly_sources.bed"),
        partners=workflow.source_path("../resources/dux4/chm13v2.0.partners.bed"),
        helper=workflow.source_path("../scripts/dux4_local_assembly.py"),
    output:
        reads=DUX4_ASM_DIR + "/target_reads.fastq.gz",
        stats=DUX4_ASM_DIR + "/extraction_stats.tsv",
    params:
        padding=DUX4_ASM.get("partner_padding", 10000),
    threads: 4
    resources:
        mem_mb=4096,
        disk_mb=8192,
    log:
        "logs/{sample}/dux4_local_extract.log",
    benchmark:
        "benchmarks/{sample}/dux4_local_extract.tsv"
    conda:
        "../envs/dux4_assembly_io.yaml"
    shell:
        """
        python {input.helper:q} extract --bam {input.bam:q} --index {input.index:q} \
            --fai {input.fai:q} --sources {input.sources:q} --partners {input.partners:q} \
            --partner-padding {params.padding} --threads {threads} \
            --fastq {output.reads:q} --stats {output.stats:q} > {log:q} 2>&1
        """


rule local_assemble_dux4:
    input:
        reads=rules.extract_dux4_igh_reads.output.reads,
        helper=workflow.source_path("../scripts/dux4_local_assembly.py"),
    output:
        contigs=DUX4_ASM_DIR + "/contigs.fa",
        status=DUX4_ASM_DIR + "/assembly_status.tsv",
    params:
        ont_mode=DUX4_ASM.get("ont_mode", "nano-hq"),
        min_reads=DUX4_ASM.get("min_reads", 3),
        min_read_length=DUX4_ASM.get("min_read_length", 1000),
    threads: 12
    resources:
        mem_mb=16000,
        disk_mb=16000,
    log:
        "logs/{sample}/dux4_local_assemble.log",
    benchmark:
        "benchmarks/{sample}/dux4_local_assemble.tsv"
    conda:
        "../envs/dux4_assembly_ont.yaml"
    shell:
        """
        python {input.helper:q} assemble --fastq {input.reads:q} \
            --platform ont --ont-mode {params.ont_mode:q} --threads {threads} \
            --min-reads {params.min_reads} --min-read-length {params.min_read_length} \
            --contigs {output.contigs:q} --status {output.status:q} > {log:q} 2>&1
        """


rule align_contigs_to_ref:
    input:
        contigs=rules.local_assemble_dux4.output.contigs,
        ref=config["fasta"],
        fai=rules.faidx.output.fai,
        sources=workflow.source_path("../resources/dux4/chm13v2.0.assembly_sources.bed"),
        partners=workflow.source_path("../resources/dux4/chm13v2.0.partners.bed"),
        helper=workflow.source_path("../scripts/dux4_local_assembly.py"),
    output:
        bam=DUX4_ASM_DIR + "/contigs_aligned.bam",
        bai=DUX4_ASM_DIR + "/contigs_aligned.bam.bai",
        summary=DUX4_ASM_DIR + "/breakpoints.tsv",
        bed=DUX4_ASM_DIR + "/breakpoints.bed",
    params:
        padding=DUX4_ASM.get("partner_padding", 10000),
        anchor=DUX4_ASM.get("min_anchor_bases", 100),
        mapq=DUX4_ASM.get("min_partner_mapq", 20),
        gap=DUX4_ASM.get("max_query_gap", 1000),
        overlap=DUX4_ASM.get("max_query_overlap", 100),
    threads: 8
    resources:
        mem_mb=16000,
        disk_mb=4096,
    log:
        "logs/{sample}/dux4_local_align.log",
    benchmark:
        "benchmarks/{sample}/dux4_local_align.tsv"
    conda:
        "../envs/dux4_assembly_io.yaml"
    shell:
        """
        python {input.helper:q} align --contigs {input.contigs:q} --ref {input.ref:q} \
            --fai {input.fai:q} --sources {input.sources:q} --partners {input.partners:q} \
            --partner-padding {params.padding} --threads {threads} --preset asm20 \
            --min-anchor-bases {params.anchor} --min-partner-mapq {params.mapq} \
            --max-query-gap {params.gap} --max-query-overlap {params.overlap} \
            --bam {output.bam:q} --index {output.bai:q} \
            --summary {output.summary:q} --bed {output.bed:q} > {log:q} 2>&1
        """
