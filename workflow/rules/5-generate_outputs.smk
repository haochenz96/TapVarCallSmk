rule bcftools_pass2_filter_sc_f_m2_calls:
    # for each single cell VCF, hard filter (v3) and index
    input:
        expand("mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz", sample_name = sample_name, cell_barcode = bars),
    output:
        expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
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
            {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

            tabix {params.bcf_outpath}/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

        done
        """

rule bcftools_pass2_merge_m2_sc_f_calls:
    # merge force-called single-cell VCFs; index; generate statistics.
    input:
        filtered_sc_f_vcfs = expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz", cell_barcode = bars, sample_name = sample_name),
        filtered_sc_f_vcfs_index = expand("bcf_pass2/filtered_vcfs/{sample_name}_{cell_barcode}_hard_filtered.vcf.gz.tbi", cell_barcode = bars, sample_name = sample_name),
    output:
        merged_prev_filtered_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF_AF_DP_filters.prev_filtered.vcf.gz", sample_name = sample_name),
        merged_prev_filtered_vcf_stats = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF_AF_DP_filters.prev_filtered.vcf.gz.stats", sample_name = sample_name),
    conda:
        "../envs/bcftools.yaml"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools merge --merge none {input.filtered_sc_f_vcfs} | \
        bcftools sort | \
        bcftools filter -i 'N_PASS(GT!="mis") > 2' | \
        bgzip -c > \
        {output.merged_prev_filtered_vcf}

        tabix {output.merged_prev_filtered_vcf}

        bcftools stats {output.merged_prev_filtered_vcf} > {output.merged_prev_filtered_vcf_stats}
        """

# rule write_h5:
#     # merge force-called single-cell VCFs; index; generate statistics.
#     input:
#         merged_f_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_f_VCF.vcf.gz", sample_name = sample_name),
#         read_counts_tsv = expand("tap_pipeline_output")
#     output:
#         output_h5 = expand("OUTPUTS/{sample_name}_DNA_CNV.h5", sample_name = sample_name)
#     params:
#         metadata = {
#             'sample_name': sample_name,
#             'genome_version': config['genome_version'],
#             'panel_version': config['panel_version'],
#             'filterm2_ops': config['mutect2']['filterm2_ops'],
#             'bcftools_filters': config['bcftools']['filters']
#         },
#         metadata_json = json.dumps(metatdata)
#     conda:
#         "../envs/mosaic.yaml"
#     threads: 4
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 2000,
#         time_min = 59
#     script:
#         "../scripts/write_h5.py"
