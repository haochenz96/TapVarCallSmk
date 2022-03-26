rule bcftools_pass2_normalize_sc_f_m2_calls:
    # for each single cell VCF, apply only DP filter; then split multiallelic sites and normalize so that they can be merged in the next step
    input:
        expand("mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz", sample_name = sample_name, cell_barcode = bars),
    output:
        expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
    params:
        bcf_filters = parse_bcf_filters(config['bcftools']['pass2_filters']),
        bcf_outpath = "bcf_pass2/filtered_vcfs"
    log:
        "logs/bcftools_pass2/filter_sc_f_m2_calls.log"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        for m2_vcf in {input}; do 
            
            out_name=$(basename $m2_vcf)
            echo 'processing ${{out_name}}'
            
            bcftools norm -m- $m2_vcf | \
            bcftools filter -i '{params.bcf_filters}' | \
            bgzip -c > \
            {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/filtered.vcf.gz}}

            tabix {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/filtered.vcf.gz}}

        done
        """

rule bcftools_pass2_merge_m2_sc_f_calls:
    # merge force-called single-cell VCFs; index; generate statistics.
    input:
        filtered_sc_f_vcfs = expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        filtered_sc_f_vcfs_index = expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
    output:
        merged_filtered_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz", sample_name = sample_name),
        merged_filtered_vcf_stats = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz.stats", sample_name = sample_name),
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools merge --merge none {input.filtered_sc_f_vcfs} | \
        bcftools sort | \
        bgzip -c > \
        {output.merged_filtered_vcf}

        tabix {output.merged_filtered_vcf}

        bcftools stats {output.merged_filtered_vcf} > {output.merged_filtered_vcf_stats}
        """

rule bcftools_pass2_intersect_f_q_vcfs:
    # intersect q_vcf and f_vcf
    input:
        f_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz", sample_name = sample_name),
        q_vcf = expand("bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz", sample_name = sample_name),
    output:
        output_vcf = expand("bcf_pass2/{sample_name}-f_q_intersected.vcf.gz", sample_name = sample_name),
        output_vcf_stats = expand("bcf_pass2/{sample_name}-f_q_intersected.vcf.gz.stats", sample_name = sample_name),
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools isec \
            -n=2 \
            -w2 \
            -Oz \
            -o {output.output_vcf} \
            {input.q_vcf} \
            {input.f_vcf}
        
        tabix {output.output_vcf}
        bcftools stats {output.output_vcf} > {output.output_vcf_stats}
        """