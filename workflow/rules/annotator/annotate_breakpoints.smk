rule annotate_breakpoints:
    input:
        vcf=rules.survivor_final.output.vcf,
        genes=config["snpeff_build"]["genes"],
    output:
        "survivor/{sample}/final/{sample}.{type_sv}.breakpoints.tsv",
    script:
        "../../scripts/annotate_breakpoints.py"
