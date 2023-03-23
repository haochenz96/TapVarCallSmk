rule step6_annotate_h5_with_bulk:
    # annotate H5 file
    # ---- get CRAVAT annotation to SNV's of interest
    input:
        input_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.h5",
        bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.h5",
        ponv1_filtered_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5",
        tlod_metrics_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.tsv",
        tlod_metrics_ponv1_filtered_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.ponv1_filtered.tsv",
    params:
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log"
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP6-annotate_h5_with_bulk.py"

rule step6_query_h5_with_CRAVAT:
    input:
        ponv1_filtered_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5",
        # bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        cravat_output = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.txt",
    params:
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log",
        opencravat_username = config['CRAVAT']['username'],
        opencravat_password = config['CRAVAT']['password'],
        cravat_annotators = config['CRAVAT']['annotators'],
        cravat_input_path = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_input.txt",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP6-query_h5_with_CRAVAT.py"

rule step6_annotate_h5_with_CRAVAT:
    input:
        ponv1_filtered_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5",
        cravat_output = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.txt",
    output:
        cravat_annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.cravat_annotated.h5",
        cravat_output_cleaned = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.cleaned.txt",
    params:
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP6-annotate_h5_with_CRAVAT.py"