rule restore_annotsv_length:
    input:
        annotsv="{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv",
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
    output:
        tsv="{caller}/{sample}/{caller}.{type_sv}.annotsv.length.tsv",
    script:
        "../../scripts/restore_annotsv_length.py"
