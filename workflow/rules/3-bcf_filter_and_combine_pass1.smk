rule bcftools_pass_1_filter_sc_m2_calls:
    # for each single cell VCF, hard filter (v3) and index
    input:
        expand("mutect2_sc_pass1/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz", cell_barcode = bars, sample_name = sample_name),
    output:
        expand("bcf_pass1/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        expand("bcf_pass1/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
    params:
        bcf_filters = parse_bcf_filters(config['bcftools']['pass1_filters']),
        bcf_outpath = "bcf_pass1/filtered_vcfs"
    log:
        "logs/bcftools/filter_sc_m2_calls.log"
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
            {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

            tabix {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

        done
        """
    
rule bcftools_pass1_merge_filtered_m2_calls:
    # merge filtered single cell VCFs; filter again based on mutational prevalence across all cells; create statistics.
    input:
        filtered_sc_vcfs = expand("bcf_pass1/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        filtered_sc_vcfs_index = expand("bcf_pass1/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
    output:
        merged_prev_filtered_vcf = expand("bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz", sample_name = sample_name),
        merged_prev_filtered_vcf_stats = expand("bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz.stats", sample_name = sample_name),
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools merge --merge none {input.filtered_sc_vcfs} | \
        bcftools sort | \
        bcftools filter -i 'N_PASS(GT!="mis") > 2' | \
        bgzip -c > \
        {output.merged_prev_filtered_vcf}

        tabix {output.merged_prev_filtered_vcf}

        bcftools stats {output.merged_prev_filtered_vcf} > {output.merged_prev_filtered_vcf_stats}
        """
    




