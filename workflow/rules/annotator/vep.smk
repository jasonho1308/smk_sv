rule vep:
    conda:
        "../../envs/vep.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        cache=path_cache_vep,
    output:
        vcf=protected("{caller}/{sample}/{caller}.{type_sv}.vep.vcf"),
        html="{caller}/{sample}/{caller}.{type_sv}.vep.html",
    params:
        cache=lambda wildcards, input: str(Path(input.cache).parents[1]),
        version=config["version_vep"],
        genome=config["genome"],
        species=config["species"],
        max_size=config["max_size_vep"],
        arg_buffer_size=lambda wildcards: (
            "--buffer_size 50" if wildcards.type_sv == "BND" else ""
        ),
    threads: 1
    log:
        "logs/{sample}/vep.{caller}.{type_sv}.log",
    shell:
        """
        # VEP reads the original input VCF directly; do not create unused intermediate VCFs.
        vep \\
            -i {input.vcf} -o {output.vcf} --stats_file {output.html} \\
            --species {params.species} --assembly {params.genome} \\
            --cache_version {params.version} --fasta {input.fasta} \\
            --dir_cache {params.cache} --fork {threads} \\
            --force_overwrite --cache --vcf --everything --filter_common \\
            --per_gene --total_length --offline --format vcf --dont_skip \\
            --max_sv_size {params.max_size} \\
            {params.arg_buffer_size} \\
            1> {log} 2>&1
        """
