# *--------------------------------------------------------------------------* #
# * Configuration                                                            * #
# *--------------------------------------------------------------------------* #
from snakemake.utils import validate


include: "utils.smk"
include: "constants.smk"


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.json")


# Keep the rule interface compact while exposing a grouped user configuration.
config = {
    "dir_run": config["run"]["dir"],
    "mapper": config["workflow"]["mapper"],
    "callers": config["workflow"]["callers"],
    "annotators": config["workflow"]["annotators"],
    **config["inputs"],
    **config["filtering"],
    "terms_relative": config["annotation"]["terms_relative"],
    "species": config["annotation"]["species"],
    "genome": config["annotation"]["genome"],
    "version_snpeff": config["annotation"]["snpeff"]["version"],
    "cache_snpeff": config["annotation"]["snpeff"]["cache"],
    "version_vep": config["annotation"]["vep"]["version"],
    "cache_vep": config["annotation"]["vep"]["cache"],
    "max_size_vep": config["annotation"]["vep"]["max_size"],
    "version_annotsv": config["annotation"]["annotsv"]["version"],
    "cache_annotsv": config["annotation"]["annotsv"]["cache"],
    "bed_tandem_repeats": config["caller_settings"]["tandem_repeats"],
    "bed_nvtr": config["caller_settings"]["nvtr"],
    "config_nanosv": config["caller_settings"]["nanosv"]["config"],
    "bed_nanosv": config["caller_settings"]["nanosv"].get("bed"),
    "model_clair3": config["caller_settings"]["clair3"]["model"],
    "model_svision": config["caller_settings"]["svision"]["model"],
    "resource_downloads": config["resources"]["downloads"],
}


pepfile: "config/pep/config.yaml"


if config["dir_run"] and config["dir_run"] is not None:

    workdir: config["dir_run"]


# *--------------------------------------------------------------------------* #
# * Constant-like variables                                                  * #
# *--------------------------------------------------------------------------* #
SAMPLES = pep.sample_table["sample_name"]
SPECIES = config["species"]
CALLERS = sorted(config["callers"])
MAPPER = config["mapper"]
ANNOTATORS = config["annotators"]


# *--------------------------------------------------------------------------* #
# * Wildcard constraints                                                     * #
# *--------------------------------------------------------------------------* #
wildcard_constraints:
    sample=r"|".join(SAMPLES),
    type_sv=r"|".join(TYPES_SV),
    caller=r"|".join(CALLERS),


# *--------------------------------------------------------------------------* #
# * Files and directories required by rules                                  * #
# *--------------------------------------------------------------------------* #
path_cache_snpeff = (
    f"{config['cache_snpeff']}/{config['genome']}.{config['version_snpeff']}"
)
path_cache_vep = f"{config['cache_vep']}/{config['species']}/{config['version_vep']}_{config['genome']}"
svision_model_files = [
    f"{config['model_svision']}.data-00000-of-00001",
    f"{config['model_svision']}.index",
    f"{config['model_svision']}.meta",
]

vcfs_svision = multiext(
    "svision/{sample}/chroms/{sample}",
    *[f".{chrom}.svision.s{config['min_reads']}.graph.vcf" for chrom in CHROMS],
)
