from shlex import quote


rule vep:
    conda:
        "../../envs/vep.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        fai=f"{config['fasta']}.fai",
        **(
            {"cache": path_cache_vep}
            if config["vep_mode"] == "cache"
            else {"annotation": path_vep_annotation, "index": f"{path_vep_annotation}.tbi"}
        ),
    output:
        vcf=protected("{caller}/{sample}/{caller}.{type_sv}.vep.vcf"),
        html="{caller}/{sample}/{caller}.{type_sv}.vep.html",
    params:
        source_args=(
            f"--offline --cache --cache_version {config['version_vep']} "
            f"--dir_cache {quote(config['cache_vep'])} --everything --filter_common"
            if config["vep_mode"] == "cache"
            else f"--{config['vep_mode']} {quote(path_vep_annotation)}"
        ),
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
        vep \\
            -i {input.vcf:q} -o {output.vcf:q} --stats_file {output.html:q} \\
            --species {params.species:q} --assembly {params.genome:q} \\
            --fasta {input.fasta:q} --fork {threads} \\
            {params.source_args} \\
            --force_overwrite --vcf --format vcf --dont_skip \\
            --symbol --biotype --canonical --numbers --hgvs --protein \\
            --transcript_version --variant_class --per_gene --total_length \\
            --max_sv_size {params.max_size} \\
            {params.arg_buffer_size} \\
            > {log:q} 2>&1
        """
