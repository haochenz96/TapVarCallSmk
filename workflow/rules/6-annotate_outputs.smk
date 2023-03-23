rule step6_annotate_h5_with_bulk:
    # annotate H5 file
    # ---- get CRAVAT annotation to SNV's of interest
    input:
        input_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.h5",
        bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        bulk_annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.h5",
        # ponv1_filtered_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5",
        # tlod_metrics_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.tsv",
        # tlod_metrics_ponv1_filtered_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.ponv1_filtered.tsv",
    params:
        bulk_annotate_script = '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/scripts/STEP6-annotate_h5_with_bulk.py'
    #     log_file = "{sample_name}/logs/{sample_name}.snakemake.STEP6.log"
    log:
        std = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}.bulk_annotation.log",
        err = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}.bulk_annotation.err",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt*4000,
        time_min = 59
    shell:
        # "../scripts/STEP6-annotate_h5_with_bulk.py"
        """
        python {params.bulk_annotate_script} \
            --sample_name {wildcards.sample_name} \
            --input_h5 {input.input_h5} \
            --bulk_info_yaml {input.bulk_info_yaml} \
            --output_bulk_annotated_h5 {output.bulk_annotated_h5} \
            --output_ponv1_filtered_h5 {output.ponv1_filtered_h5} \
            1> {log.std} 2> {log.err}
        """

rule step6_query_h5_with_CRAVAT:
    input:
        ponv1_filtered_h5 = ancient("{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5"),
        # bulk_info_yaml = "{sample_name}/reference/{sample_name}_matched_bulk_info.yaml",
    output:
        cravat_output = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.txt",
    params:
        cravat_query_script = '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/scripts/STEP6-query_h5_with_CRAVAT.py',
        cravat_settings_yaml = "/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/configs/cravat_settings.yaml",
        output_dir = "{sample_name}/OUTPUTS_from_m2_f/annotations",
    conda:
        "../envs/mosaic-custom.yaml"
    log:
        std = "{sample_name}/logs/{sample_name}.snakemake.cravat.log",
        err = "{sample_name}/logs/{sample_name}.snakemake.cravat.err",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        python {params.cravat_query_script} \
            --sample_name {wildcards.sample_name} \
            --input_h5 {input.ponv1_filtered_h5} \
            --cravat_settings_yaml {params.cravat_settings_yaml} \
            --output_dir {params.output_dir} \
            1> {log.std} 2> {log.err}
        """

rule step6_annotate_h5_with_CRAVAT:
    input:
        # every time the h5 is read, it gets a new timestamp and therefore needs to be handled by ancient()
        ponv1_filtered_h5 = ancient("{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5"),
        cravat_output = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.txt",
    output:
        cravat_annotated_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}.cravat_annotated.h5",
        cravat_output_cleaned = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.cleaned.txt",
    params:
        cravat_annotate_script = '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/scripts/STEP6-annotate_h5_with_CRAVAT.py',
        output_dir = "{sample_name}/OUTPUTS_from_m2_f",
    log:
        std = "{sample_name}/logs/{sample_name}.snakemake.cravat.log",
        err = "{sample_name}/logs/{sample_name}.snakemake.cravat.err",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: max(10 * input.size_mb * attempt, 4000),
        time_min = 59
    shell:
        """
        python {params.cravat_annotate_script} \
            --sample_name {wildcards.sample_name} \
            --input_h5 {input.ponv1_filtered_h5} \
            --cravat_report {input.cravat_output} \
            --output_dir {params.output_dir} \
            1> {log.std} 2> {log.err}
        """