rule bcftools_merge_m2_sc_f_calls:
    # merge force-called single-cell VCFs; index; generate statistics.
    input:
        expand("mutect2_sc_f_pass2/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz", sample_name = sample_name, cell_barcode = bars),
    output:
        merged_f_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_f_VCF.vcf.gz", sample_name = sample_name),
        merged_f_vcf_stats = expand("bcf_pass2/combined_vcf/{sample_name}-combined_f_VCF.vcf.gz.stats", sample_name = sample_name),
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools merge --merge none {input} | \
        bcftools sort | \
        bgzip -c > \
        {output.merged_f_vcf}

        tabix {output.merged_f_vcf}

        bcftools stats {output.merged_f_vcf} > {output.merged_f_vcf_stats}
        """

rule bcftools_merge_m2_sc_f_calls:
    # merge force-called single-cell VCFs; index; generate statistics.
    input:
        merged_f_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_f_VCF.vcf.gz", sample_name = sample_name),
    output:
        output_h5 = expand("OUTPUTS/{sample_name}_DNA_CNV.h5", sample_name = sample_name)
    conda:
        "../envs/mosaic.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/write_h5.py"
