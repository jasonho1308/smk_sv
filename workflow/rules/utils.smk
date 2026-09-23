def get_resource_config(name):
    assembly = {"hg19": "GRCh37", "hg38": "GRCh38"}.get(
        config["genome"], config["genome"]
    )
    downloads = config.get("resource_downloads", {})
    resource = downloads.get(name)
    if resource is None and name == "nvtr":
        resource = downloads.get("tandem_repeats", {})
    resource = resource or {}

    if "assemblies" in resource:
        return resource["assemblies"].get(assembly, {})
    return resource


def get_targets(wildcards=None):
    targets = []

    targets += [
        f"survivor/{sample}/final/{sample}.{type_sv}.merged.vcf"
        for sample in SAMPLES
        for type_sv in TYPES_SV
    ]

    return targets


def get_annotsv_cache_outputs():
    if SPECIES in ["homo_sapiens"]:
        exomiser_version = get_resource_config("annotsv_exomiser")["version"]
        return {
            "human": f"{config['cache_annotsv']}/Annotations_Human",
            "exomiser": f"{config['cache_annotsv']}/Annotations_Exomiser/{exomiser_version}",
        }
    else:
        raise ValueError("Unsupported species")


def get_resource_targets(wildcards=None):
    "Define the targets for the resource.smk."
    targets = [config["fasta"], f"{config['fasta']}.fai", config["index_minimap2"]]

    if "sniffles" in CALLERS:
        targets.append(config["bed_tandem_repeats"])
    if "severus" in CALLERS:
        targets.extend([config["bed_nvtr"], config["model_clair3"]])
    if "clair3" in CALLERS:
        targets.append(config["model_clair3"])
    if "svision" in CALLERS:
        targets.extend(svision_model_files)
    if "nanosv" in CALLERS:
        targets.append(config["config_nanosv"])
        if config["bed_nanosv"]:
            targets.append(config["bed_nanosv"])
    if "snpeff" in ANNOTATORS:
        targets.append(path_cache_snpeff)
    if "vep" in ANNOTATORS:
        targets.append(path_cache_vep)
    if "annotsv" in ANNOTATORS:
        targets.extend(get_annotsv_cache_outputs().values())

    return list(dict.fromkeys(targets))


def get_convert_snpeff_arguments(wildcards):
    caller = wildcards.caller

    fields = CALLER2FMTS.get(caller)
    if fields is None:
        raise ValueError("Unsupported caller")

    arg = " ".join(f"GEN[*].{field}" for field in fields)

    return arg
