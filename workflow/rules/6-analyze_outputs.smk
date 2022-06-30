rule step6_annotate_h5_with_bulk:
    # annotate H5 file
    # ---- get CRAVAT annotation to SNV's of interest
    input:
        input_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.h5",
        bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.h5",
        tlod_metrics_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.tsv",
    params:
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log"
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP6-annotate_vars_with_bulk.py"

rule step6_annotate_h5_with_CRAVAT:
    input:
        bulk_annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.h5",
        # bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        cravat_df = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT.txt",
    params:
        cravat_annotators = config['CRAVAT']['annotators'],
        log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log"
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    script:
        "../scripts/STEP6-annotate_h5_with_CRAVAT.py"