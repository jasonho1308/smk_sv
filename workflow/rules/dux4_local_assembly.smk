# Optional local assembly, independent of the read-level screen and SV consensus.
from snakemake.exceptions import WorkflowError

DUX4_ASM = config["dux4_local_assembly"]
DUX4_ASM_DIR = "results/local_asm/{sample}"
DUX4_ASM_REF = DUX4_ASM.get("ref", config["fasta"])
DUX4_ASM_PLATFORM = DUX4_ASM.get("platform", "ont")
DUX4_ASM_CONTAINER = DUX4_ASM.get("container", "docker://condaforge/miniforge3:25.3.1-0")
DUX4_ASM_REGIONS = DUX4_ASM.get("regions", {})


def dux4_asm_fai(wildcards):
    if DUX4_ASM_REF == config["fasta"]:
        return rules.faidx.output.fai
    return rules.index_dux4_local_reference.output.fai


def dux4_asm_region(kind):
    def resolve(wildcards):
        if kind in DUX4_ASM_REGIONS:
            return DUX4_ASM_REGIONS[kind]
        if config["genome"] != "CHM13" or DUX4_ASM_REF != config["fasta"]:
            raise WorkflowError("Set caller_settings.dux4_local_assembly.regions.sources and .partners "
                                "for references other than the workflow's CHM13v2.0 FASTA")
        resource = {"sources": "chm13v2.0.assembly_sources.bed", "partners": "chm13v2.0.partners.bed"}
        return workflow.source_path("../resources/dux4/" + resource[kind])
    return resolve


def dux4_asm_bam(wildcards):
    pattern = DUX4_ASM.get("bam_pattern", rules.minimap2.output.bam)
    return pattern.replace("{sample}", wildcards.sample)


def dux4_asm_bam_index(wildcards):
    if "bam_pattern" not in DUX4_ASM:
        return rules.minimap2.output.bai.replace("{sample}", wildcards.sample)
    bam = dux4_asm_bam(wildcards)
    return DUX4_ASM.get("bam_index_pattern", bam + ".bai").replace("{sample}", wildcards.sample)


def dux4_asm_targets():
    return expand(DUX4_ASM_DIR + "/{product}", sample=SAMPLES,
                  product=["assembly_status.tsv", "contigs_aligned.bam", "contigs_aligned.bam.bai",
                           "breakpoints.tsv", "breakpoints.bed"])


rule dux4_local_assembly:
    input:
        dux4_asm_targets(),


rule index_dux4_local_reference:
    input:
        ref=DUX4_ASM_REF,
    output:
        fai="results/local_asm/reference.fai",
    threads: 1
    resources:
        mem_mb=512,
        disk_mb=256,
    log:
        "logs/dux4_local_reference.log",
    benchmark:
        "benchmarks/dux4_local_reference.tsv"
    conda:
        "../envs/dux4_assembly_io.yaml"
    container:
        DUX4_ASM_CONTAINER
    shell:
        "samtools faidx --fai-idx {output.fai:q} {input.ref:q} > {log:q} 2>&1"


rule extract_dux4_igh_reads:
    # Historical rule name: collects reads at all configured partners, not just IGH.
    input:
        bam=dux4_asm_bam,
        index=dux4_asm_bam_index,
        fai=dux4_asm_fai,
        sources=dux4_asm_region("sources"),
        partners=dux4_asm_region("partners"),
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
    container:
        DUX4_ASM_CONTAINER
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
        platform=DUX4_ASM_PLATFORM,
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
        "../envs/dux4_assembly_" + DUX4_ASM_PLATFORM + ".yaml"
    container:
        DUX4_ASM_CONTAINER
    shell:
        """
        python {input.helper:q} assemble --fastq {input.reads:q} \
            --platform {params.platform:q} --ont-mode {params.ont_mode:q} --threads {threads} \
            --min-reads {params.min_reads} --min-read-length {params.min_read_length} \
            --contigs {output.contigs:q} --status {output.status:q} > {log:q} 2>&1
        """


rule align_contigs_to_ref:
    input:
        contigs=rules.local_assemble_dux4.output.contigs,
        ref=DUX4_ASM_REF,
        fai=dux4_asm_fai,
        sources=dux4_asm_region("sources"),
        partners=dux4_asm_region("partners"),
        helper=workflow.source_path("../scripts/dux4_local_assembly.py"),
    output:
        bam=DUX4_ASM_DIR + "/contigs_aligned.bam",
        bai=DUX4_ASM_DIR + "/contigs_aligned.bam.bai",
        summary=DUX4_ASM_DIR + "/breakpoints.tsv",
        bed=DUX4_ASM_DIR + "/breakpoints.bed",
    params:
        preset="asm5" if DUX4_ASM_PLATFORM == "hifi" else "asm20",
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
    container:
        DUX4_ASM_CONTAINER
    shell:
        """
        python {input.helper:q} align --contigs {input.contigs:q} --ref {input.ref:q} \
            --fai {input.fai:q} --sources {input.sources:q} --partners {input.partners:q} \
            --partner-padding {params.padding} --threads {threads} --preset {params.preset:q} \
            --min-anchor-bases {params.anchor} --min-partner-mapq {params.mapq} \
            --max-query-gap {params.gap} --max-query-overlap {params.overlap} \
            --bam {output.bam:q} --index {output.bai:q} \
            --summary {output.summary:q} --bed {output.bed:q} > {log:q} 2>&1
        """
