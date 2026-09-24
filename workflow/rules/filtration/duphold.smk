rule duphold:
    conda:
        "../../envs/duphold.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.vcf",
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        bcf=temp("{caller}/{sample}/{sample}.sorted.bcf"),
        vcf="{caller}/{sample}/{caller}.duphold.vcf",
    threads: 1
    log:
        "logs/{sample}/duphold.{caller}.log",
    shell:
        """
        {{ # Discard malformed POS=0 records before BCF conversion (duphold crashes on them).
        awk -F '\t' '/^#/ {{print; next}} $2 !~ /^[0-9]+$/ || $2+0 < 1 {{invalid++; next}} {{print}} END {{if (invalid) print "Skipped " invalid " VCF records with invalid POS" > "/dev/stderr"}}' {input.vcf} \\
            | bcftools sort -Ou - > {output.bcf}

        export DUPHOLD_SAMPLE_NAME={wildcards.sample}
        duphold \\
            -t {threads} \\
            -v {output.bcf} \\
            -b {input.bam} \\
            -f {input.fasta} \\
            -o {output.vcf}; }} \\
        1> {log} 2>&1
        """
