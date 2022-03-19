rule write_h5:
    # merge force-called single-cell VCFs; index; generate statistics.
    input:
        merged_prev_filtered_vcf = expand("bcf_pass2/combined_vcf/{sample_name}-combined_VCF_DP_AF_filters.prev_filtered.vcf.gz", sample_name = sample_name),
        read_counts_tsv = expand("tap_pipeline_output/results/tsv/{sample_name}_{genome_version}.tube1.barcode.cell.distribution.tsv", sample_name = sample_name, genome_version = genome_version)
    output:
        output_h5 = expand("OUTPUTS/{sample_name}_DNA_CNV.h5", sample_name = sample_name)
    params:
        metadata_json = json.dumps({
            "sample_name": sample_name,
            "genome_version": config['reference_info']['genome_version'],
            "panel_version": config['reference_info']['panel_version'],
            "filterm2_ops": config['mutect2']['filterm2_ops'],
            "bcftools_filters": config['bcftools']['pass2_filters']
        }),
    conda:
        "../envs/mosaic.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/write_h5.py"
