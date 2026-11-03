# Custom builds and prebuilt downloads are mutually exclusive cache producers.
if config["snpeff_build"]:

    rule download_snpeff_genes:
        conda:
            "../envs/resources.yaml"
        output:
            protected(config["snpeff_build"]["genes"]),
        params:
            url=get_resource_config("snpeff_genes").get("url", ""),
            checksum=get_resource_config("snpeff_genes").get("checksum", ""),
            compression=get_resource_config("snpeff_genes").get("compression", "none"),
            postprocess="none",
        log:
            "logs/resources/snpeff_genes.log",
        script:
            "../scripts/download_resource.py"


    rule build_snpeff_cache:
        conda:
            "../envs/snpeff.yaml"
        input:
            fasta=config["fasta"],
            fai=f"{config['fasta']}.fai",
            genes=config["snpeff_build"]["genes"],
            cds=[config["snpeff_build"]["cds"]] if config["snpeff_build"].get("cds") else [],
            protein=[config["snpeff_build"]["protein"]] if config["snpeff_build"].get("protein") else [],
        output:
            cache=protected(directory(path_cache_snpeff)),
        params:
            database=f"{config['genome']}.{config['version_snpeff']}",
            format=config["snpeff_build"].get("format", "gff3"),
        resources:
            mem_mb=32000,
        log:
            "logs/resources/snpeff_build.log",
        script:
            "../scripts/build_snpeff_database.py"

else:

    rule download_snpeff_cache:
        conda:
            "../envs/awscli.yaml"
        output:
            cache=protected(directory(path_cache_snpeff)),
        params:
            source=f"{config['genome']}.{config['version_snpeff']}",
            destination=lambda wildcards, output: str(Path(output.cache).parent),
            cache=lambda wildcards, output: str(Path(output.cache).parent),
        log:
            "logs/resources/snpeff_cache.log",
        shell:
            r"""
            set -euo pipefail
            if [[ {params.source:q} == CHM13.* ]]; then
                echo "No official CHM13 SnpEff database is available. Configure annotation.snpeff.build; see README.md." >&2
                exit 1
            fi
            url_base="s3://annotation-cache/snpeff_cache"
            mkdir -p {params.destination:q}
            if aws s3 ls --no-sign-request "$url_base/" \
                | awk '$1 == "PRE" {{sub(/\/$/, "", $2); print $2}}' \
                | grep -Fxq {params.source:q}; then
                aws s3 --no-sign-request sync \
                    "$url_base/{params.source}" {params.destination:q}
            else
                snpEff download {params.source:q} -dataDir {params.cache:q}
            fi > {log:q} 2>&1
            """
