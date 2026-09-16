# Download and prepare versioned reference data used by the selected workflow tools.

localrules: resources


rule resources:
    input:
        get_resource_targets,


rule download_reference:
    conda:
        "../envs/resources.yaml"
    output:
        protected(config["fasta"]),
    params:
        url=get_resource_config("reference").get("url", ""),
        checksum=get_resource_config("reference").get("checksum", ""),
        compression=get_resource_config("reference").get("compression", "gzip"),
        postprocess=get_resource_config("reference").get("postprocess", "none"),
    log:
        "logs/resources/reference.log",
    script:
        "../scripts/download_resource.py"


rule download_tandem_repeats:
    conda:
        "../envs/resources.yaml"
    output:
        protected(config["bed_tandem_repeats"]),
    params:
        url=get_resource_config("tandem_repeats").get("url", ""),
        checksum=get_resource_config("tandem_repeats").get("checksum", ""),
        compression=get_resource_config("tandem_repeats").get("compression", "none"),
        postprocess=get_resource_config("tandem_repeats").get("postprocess", "none"),
    log:
        "logs/resources/tandem_repeats.log",
    script:
        "../scripts/download_resource.py"


if config["bed_nvtr"] != config["bed_tandem_repeats"]:

    rule download_nvtr:
        conda:
            "../envs/resources.yaml"
        output:
            protected(config["bed_nvtr"]),
        params:
            url=get_resource_config("nvtr").get("url", ""),
            checksum=get_resource_config("nvtr").get("checksum", ""),
            compression=get_resource_config("nvtr").get("compression", "none"),
            postprocess=get_resource_config("nvtr").get("postprocess", "none"),
        log:
            "logs/resources/nvtr.log",
        script:
            "../scripts/download_resource.py"


rule download_nanosv_config:
    conda:
        "../envs/resources.yaml"
    output:
        protected(config["config_nanosv"]),
    params:
        url=get_resource_config("nanosv_config").get("url", ""),
        checksum=get_resource_config("nanosv_config").get("checksum", ""),
        compression="none",
        postprocess=get_resource_config("nanosv_config").get("postprocess", "none"),
    log:
        "logs/resources/nanosv_config.log",
    script:
        "../scripts/download_resource.py"


if config["bed_nanosv"]:
    rule download_nanosv_bed:
        conda:
            "../envs/resources.yaml"
        output:
            protected(config["bed_nanosv"]),
        params:
            url=get_resource_config("nanosv_bed").get("url", ""),
            checksum=get_resource_config("nanosv_bed").get("checksum", ""),
            compression="none",
            postprocess=get_resource_config("nanosv_bed").get("postprocess", "none"),
        log:
            "logs/resources/nanosv_bed.log",
        script:
            "../scripts/download_resource.py"


rule download_clair3_model:
    conda:
        "../envs/resources.yaml"
    output:
        model=protected(directory(config["model_clair3"])),
    params:
        url=get_resource_config("clair3_model").get("url", ""),
        checksum=get_resource_config("clair3_model").get("checksum", ""),
        archive_type="tar",
        archive_root=get_resource_config("clair3_model").get("archive_root", ""),
        destination=lambda wildcards, output: str(output.model),
    log:
        "logs/resources/clair3_model.log",
    script:
        "../scripts/install_archive.py"


rule download_svision_model:
    conda:
        "../envs/resources.yaml"
    output:
        data=protected(svision_model_files[0]),
        index=protected(svision_model_files[1]),
        meta=protected(svision_model_files[2]),
    params:
        url=get_resource_config("svision_model").get("url", ""),
        destination=lambda wildcards, output: str(Path(output.index).parent),
        prefix=lambda wildcards, output: Path(output.index).stem,
    log:
        "logs/resources/svision_model.log",
    shell:
        r"""
        set -euo pipefail
        test -n {params.url:q} || {{
            echo "No download URL configured for resource_downloads.svision_model" >&2
            exit 1
        }}
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT
        gdown --folder {params.url:q} --output "$tmpdir" > {log:q} 2>&1
        mkdir -p {params.destination:q}
        for suffix in data-00000-of-00001 index meta; do
            filename="{params.prefix}.${{suffix}}"
            source=$(find "$tmpdir" -type f -name "$filename" -print -quit)
            test -n "$source" || {{
                echo "SVision model file missing from download: $filename" >&2
                exit 1
            }}
            mv "$source" "{params.destination}/$filename"
        done
        """


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


rule download_annotsv_human:
    conda:
        "../envs/resources.yaml"
    output:
        annotations=protected(directory(get_annotsv_cache_outputs()["human"])),
    params:
        url=get_resource_config("annotsv_human").get("url", ""),
        checksum=get_resource_config("annotsv_human").get("checksum", ""),
        archive_type="tar",
        archive_root=get_resource_config("annotsv_human").get(
            "archive_root", "Annotations_Human"
        ),
        destination=lambda wildcards, output: str(output.annotations),
    log:
        "logs/resources/annotsv_human.log",
    script:
        "../scripts/install_archive.py"


rule download_annotsv_exomiser:
    conda:
        "../envs/resources.yaml"
    output:
        annotations=protected(directory(get_annotsv_cache_outputs()["exomiser"])),
    params:
        url=get_resource_config("annotsv_exomiser").get("url", ""),
        checksum=get_resource_config("annotsv_exomiser").get("checksum", ""),
        archive_type="zip",
        archive_root="",
        destination=lambda wildcards, output: str(output.annotations),
    log:
        "logs/resources/annotsv_exomiser.log",
    script:
        "../scripts/install_archive.py"
