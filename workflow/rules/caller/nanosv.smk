rule nanosv:
    conda:
        "../../envs/nanosv.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        config=config["config_nanosv"],
    output:
        vcf=temp("nanosv/{sample}/tmp.vcf"),
    threads: 1
    params:
        bed=f"-b {config['bed_nanosv']}" if config["bed_nanosv"] else "",
    log:
        "logs/{sample}/nanosv.log",
    shell:
        """
        NanoSV \\
            -t {threads} \\
            -s samtools \\
            -c {input.config} \\
            {params.bed} \\
            -o {output.vcf} \\
            {input.bam} \\
            1> {log} 2>&1
        """


rule format_nanosv:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="nanosv/{sample}/nanosv.vcf",
    output:
        tab=temp("debreak/{sample}/rename.tab"),
        vcf=protected("nanosv/{sample}/nanosv.vcf"),
    params:
        header=f"{MAPPER}/{{sample}}/{{sample}}\\t{{sample}}",
    log:
        "logs/{sample}/format_nanosv.log",
    shell:
        """
        {{ echo -e "{params.header}" > {output.tab}
        bcftools reheader -s {output.tab} {input.vcf} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
