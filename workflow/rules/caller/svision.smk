rule svision:
    container:
        "docker://jiadongxjtu/svision:latest"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
        model_data=svision_model_files[0],
        model_index=svision_model_files[1],
        model_meta=svision_model_files[2],
    output:
        # SVision may omit the VCF when there are no calls; the formatter skips empty files.
        vcf=touch(
            f"svision/{{sample}}/chroms/{{chrom}}/{{sample}}.{{chrom}}.svision.s{config['min_reads']}.graph.vcf"
        ),
    params:
        dir="svision/{sample}/chroms/{chrom}",
        model=lambda wildcards, input: str(input.model_index)[: -len(".index")],
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
        max_size=config["max_size"],
    threads: 1
    log:
        "logs/{sample}/svision/{chrom}.log",
    wildcard_constraints:
        chrom=r"|".join(CHROMS),
    shell:
        """
        SVision \\
            -t {threads} \\
            -s {params.min_reads} \\
            --min_mapq {params.min_quality_mapping} \\
            --min_sv_size {params.min_size} \\
            --max_sv_size {params.max_size} \\
            --qname \\
            --graph \\
            --min_gt_depth {params.min_reads} \\
            -o {params.dir} \\
            -b {input.bam} \\
            -m {params.model} \\
            -g {input.fasta} \\
            -n {wildcards.sample}.{wildcards.chrom} \\
            -c {wildcards.chrom} \\
            1> {log} 2>&1
        """


rule format_svision:
    conda:
        "../../envs/bcftools.yaml"
    input:
        [str(vcf) for vcf in vcfs_svision],
    output:
        tab=temp("svision/{sample}/rename.tab"),
        vcf=protected("svision/{sample}/svision.vcf"),
    log:
        "logs/{sample}/format_svision.log",
    shell:
        """
        {{ lead_vcf=""
        for vcf in {input}; do
            if [ -s "$vcf" ]; then
                lead_vcf="$vcf"
                break
            fi
        done
        if [ -z "$lead_vcf" ]; then
            echo "No non-empty SVision VCF found for {wildcards.sample}" >&2
            exit 1
        fi
        lead_chrom=$(basename "$(dirname "$lead_vcf")")
        printf '%s\\t%s\\n' "{wildcards.sample}.$lead_chrom" "{wildcards.sample}" > {output.tab}

        {{ grep '^#' "$lead_vcf"; awk '!/^#/' {input}; }} \\
            | awk '/^##INFO=<ID=GFA_L/ && !f {{print "##INFO=<ID=GFA_ID,Number=.,Type=String,Description=\\"GFA_ID\\">"; f=1}} 1' \\
            | awk 'BEGIN{{FS=OFS="\\t"}} /^#/ || $5 == "<CSV>" {{print; next}} {{split($8, a, ";"); for(i in a) {{if(a[i] ~ /^SVTYPE=/) {{split(a[i], b, "="); if(b[2] == "tDUP") b[2] = "DUP:TANDEM"; gsub("<SV>", "<"b[2]">", $5)}}}}}}1' \\
            | awk '/^##ALT/ && !f {{print "##ALT=<ID=INS,Description=\\"INS\\">\\n##ALT=<ID=INV,Description=\\"INV\\">\\n##ALT=<ID=DUP,Description=\\"DUP\\">\\n##ALT=<ID=DUP:TANDEM,Description=\\"DUP:TANDEM\\">\\n##ALT=<ID=DEL,Description=\\"DEL\\">"; f=1}} 1' \\
            | awk '/^##INFO=<ID=SUPPORT,/ {{sub("Type=String", "Type=Integer")}} 1' \\
            | bcftools reheader -s {output.tab} \\
            | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
            > {output.vcf}; }} \\
        1> {log} 2>&1
        """
