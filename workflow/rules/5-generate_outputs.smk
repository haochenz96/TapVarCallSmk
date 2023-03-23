rule step5_write_h5_from_mpileup_vcf_and_rc:
    # write SNV and CNV matrices.
    input:
        # output_vcf = "{sample_name}/bcf_pass2/{sample_name}-f_q_intersected.vcf.gz",
        # output_vcf_stats = "{sample_name}/bcf_pass2/{sample_name}-f_q_intersected.vcf.gz.stats",
        output_vcf = "{sample_name}/4-bcf_genotyping/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        output_vcf_stats = "{sample_name}/4-bcf_genotyping/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz.stats",
        # the stat file ensures that the merging step finishes
        read_counts_tsv = "{sample_name}/tap_pipeline_output/results/tsv/{sample_name}.tube1.barcode.cell.distribution.tsv",
    output:
        temp_dir = temp("{sample_name}/OUTPUTS_from_mpileup/temp"),
        read_count_tsv_renamed = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.per_amplicon_read_counts.tsv",
        output_h5 = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}_DNA_CNV_bcf_mpileup.h5",
    params:
        metadata_json = json.dumps({
            "sample_name": "{sample_name}",
            "genome_version": config['reference_info']['genome_version'],
            "panel_version": config['reference_info']['panel_version'],
            "filterm2_ops": config['mutect2']['filterm2_ops'],
            "bcftools_pass1_sc_vcf_filters": config['bcftools']['pass1_sc_vcf_filters'],
            "bcftools_pass1_merged_vcf_filters": config['bcftools']['pass1_merged_vcf_filters'],
        }),
        bars_map = lambda wildcards: sample_barcode_maps[wildcards.sample_name],
    log: 
        "{sample_name}/logs/{sample_name}.snakemake.log"
    conda:
        "../envs/mosaic.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP5-write_h5_from_combined-vcf_and_rc.py"

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
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP5.log",
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
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP5-write_h5_from_combined-vcf_and_rc.py"