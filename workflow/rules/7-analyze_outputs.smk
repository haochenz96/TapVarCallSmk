rule step7-plot_tlod_metrics:
    input:
        tlod_metrics_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.tsv",
        tlod_metrics_ponv1_filtered_tsv = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_sc_TLOD_metrics.annotated.ponv1_filtered.tsv",
    output:
    
rule step7_parse_annotate_cravat_output:
    input:
        ponv1_filtered_h5 = "{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5",
        cravat_output = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_output.txt",
    output:
        
    params:
        opencravat_username = config['CRAVAT']['username'],
        opencravat_password = config['CRAVAT']['password'],
        cravat_input_path = "{sample_name}/OUTPUTS_from_m2_f/annotations/{sample_name}_CRAVAT_input.txt",
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