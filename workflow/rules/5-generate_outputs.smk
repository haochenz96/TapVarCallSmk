rule step5_write_mpileup_metadata:
    input: 
    output:
        metadata_json_f = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}_bcf_mpileup.metadata.json",
    params:
        metadata_json = json.dumps({
            "sample_name": ["{sample_name}"],
            "genome_version": config['reference_info']['genome_version'],
            "panel_version": config['reference_info']['panel_version'],
            "filterm2_ops": config['mutect2']['filterm2_ops'],
            "bcftools_pass1_sc_vcf_filters": config['bcftools']['pass1_sc_vcf_filters'],
            "bcftools_pass1_merged_vcf_filters": config['bcftools']['pass1_merged_vcf_filters'],
            "pass2_route": "bcf_mpileup; raw AD, DP extraction",
        }),
    run:
        with open(output.metadata_json_f, 'w') as f:
            f.write(params.metadata_json)


rule step5_write_h5_from_mpileup_raw_data_and_rc:
    # write SNV and CNV matrices.
    input:
        metadata_json_f = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}_bcf_mpileup.metadata.json",
        SC_MPILEUP_AD_LAYER_MERGED = "{sample_name}/4-bcf_genotyping/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "{sample_name}/4-bcf_genotyping/{sample_name}.mpileup.DP.merged.csv",
        barcode_num_map_f =  "{sample_name}/reference/{sample_name}.barcode_map.txt",
        read_counts_tsv = "{sample_name}/tap_pipeline_output/results/tsv/{sample_name}.tube1.barcode.cell.distribution.tsv", # V1
        # read_counts_tsv = lambda wildcards: config['sample_info']['sample_rc_tsvs'][wildcards.sample_name], # V2
        panel_insert_file = config['reference_info']['panel_insert_file_UCSC'], # this need to be UCSC format ('chr1' instead of '1')
        panel_amplicons_file = config['reference_info']['panel_amplicon_file'],
    output:
        read_count_tsv_renamed = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.per_amplicon_read_counts.tsv",
        output_h5 = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.mpileup.h5",
    params:
        output_dir = "{sample_name}/OUTPUTS_from_mpileup",
        write_h5_script = config['python_scripts']['WRITE_H5_FROM_MPILEUP_RAW_DATA_AND_RC'],
    log: 
        std = "{sample_name}/logs/{sample_name}.STEP5-write_h5_from_raw_data.log",
        err = "{sample_name}/logs/{sample_name}.STEP5-write_h5_from_raw_data.err",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = 59
    shell:
        """
        python {params.write_h5_script} \
            --sample_name {wildcards.sample_name} \
            --metadata_json {input.metadata_json_f} \
            --barcode_num_map_f {input.barcode_num_map_f} \
            --AD_merged_df {input.SC_MPILEUP_AD_LAYER_MERGED} \
            --DP_merged_df {input.SC_MPILEUP_DP_LAYER_MERGED} \
            --read_counts_tsv {input.read_counts_tsv} \
            --panel_insert_file {input.panel_insert_file} \
            --panel_amplicon_file {input.panel_amplicons_file} \
            --output_dir {params.output_dir} \
            --output_h5 {output.output_h5} \
            --read_counts_tsv_renamed {output.read_count_tsv_renamed} \
            1> {log.std} 2> {log.err}
        """

rule step5_write_h5_from_m2_f_vcf_and_rc:
    # write SNV and CNV matrices.
    input:
        output_vcf = "{sample_name}/4-sc_mutect_f_call/{sample_name}-f_q_intersected.vcf.gz",
        output_vcf_stats = "{sample_name}/4-sc_mutect_f_call/{sample_name}-f_q_intersected.vcf.gz.stats",
        # the stat file ensures that the merging step finishes
        read_counts_tsv = "{sample_name}/tap_pipeline_output/results/tsv/{sample_name}.tube1.barcode.cell.distribution.tsv",
    output:
        temp_dir = temp(directory("{sample_name}/OUTPUTS_from_m2_f/temp")),
        read_count_tsv_renamed = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}.per_amplicon_read_counts.tsv",
        output_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.h5",
    params:
        log_file = "{sample_name}/logs/{sample_name}.STEP5-write_h5_from_vcf.log",
        metadata_json = json.dumps({
            "sample_name": "{sample_name}",
            "genome_version": config['reference_info']['genome_version'],
            "panel_version": config['reference_info']['panel_version'],
            "filterm2_ops": config['mutect2']['filterm2_ops'],
            "bcftools_pass1_sc_vcf_filters": config['bcftools']['pass1_sc_vcf_filters'],
            "bcftools_pass1_merged_vcf_filters": config['bcftools']['pass1_merged_vcf_filters'],
        }),
        bars_map = lambda wildcards: sample_barcode_maps[wildcards.sample_name],
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: max(2 * input.size_mb * attempt, 4000),
        time_min = 59
    script:
        "../scripts/STEP5-write_h5_from_combined-vcf_and_rc.py"