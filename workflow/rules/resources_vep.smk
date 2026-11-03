if config["vep_mode"] != "cache":

    rule download_vep_genes:
        conda:
            "../envs/resources.yaml"
        output:
            protected(config["vep_annotation"]),
        params:
            url=get_resource_config("vep_genes").get("url", ""),
            checksum=get_resource_config("vep_genes").get("checksum", ""),
            compression=get_resource_config("vep_genes").get("compression", "none"),
            postprocess="none",
        log:
            "logs/resources/vep_genes.log",
        script:
            "../scripts/download_resource.py"


    rule prepare_vep_annotation:
        conda:
            "../envs/htslib.yaml"
        input:
            genes=config["vep_annotation"],
            fai=f"{config['fasta']}.fai",
        output:
            annotation=protected(path_vep_annotation),
            index=protected(f"{path_vep_annotation}.tbi"),
        threads: 1
        resources:
            mem_mb=4000,
        log:
            "logs/resources/vep_annotation.log",
        script:
            "../scripts/prepare_vep_annotation.py"

else:

    rule download_vep_cache:
        conda:
            "../envs/awscli.yaml"
        output:
            cache=protected(directory(path_cache_vep)),
        params:
            source=f"{config['version_vep']}_{config['genome']}",
            destination=lambda wildcards, output: str(Path(output.cache).parents[1]),
            cache=lambda wildcards, output: str(Path(output.cache).parents[1]),
            species=config["species"],
            version=config["version_vep"],
            genome=config["genome"],
        log:
            "logs/resources/vep_cache.log",
        shell:
            r"""
            set -euo pipefail
            if [[ {params.genome:q} == CHM13 ]]; then
                echo "No official CHM13 VEP cache is available. Use annotation.vep.mode: gff with a matching annotation; see README.md." >&2
                exit 1
            fi
            url_base="s3://annotation-cache/vep_cache"
            mkdir -p {params.destination:q}
            if aws s3 ls --no-sign-request "$url_base/" \
                | awk '$1 == "PRE" {{sub(/\/$/, "", $2); print $2}}' \
                | grep -Fxq {params.source:q}; then
                aws s3 --no-sign-request sync \
                    "$url_base/{params.source}" {params.destination:q}
            else
                vep_install \
                    --CACHEDIR {params.cache:q} \
                    --DESTDIR {params.cache:q} \
                    --PLUGINSDIR {params.cache:q}/Plugins/ \
                    --CACHE_VERSION {params.version:q} \
                    --SPECIES {params.species:q} \
                    --ASSEMBLY {params.genome:q} \
                    --PREFER_BIN --NO_UPDATE --AUTO cf
            fi > {log:q} 2>&1
            """
