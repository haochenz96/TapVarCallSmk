rule STEP7_plot_2D_TLOD_metrics:
    input:
        ponv1_filtered_h5 = ancient("{sample_name}/OUTPUTS_from_m2_f/{sample_name}_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5"),
    output:
        fig1 = "{sample_name}/OUTPUTS_from_m2_f/analysis/{sample_name}_sc-max-TLOD_scatter_2D_bulk_val.png",
        fig2 = "{sample_name}/OUTPUTS_from_m2_f/analysis/{sample_name}_sc-max-TLOD_scatter_2D_ponv2.png",
        fig3 = "{sample_name}/OUTPUTS_from_m2_f/analysis/{sample_name}_sc-max-TLOD_scatter_2D_ponv1.png",
    params:
        output_dir = "{sample_name}/OUTPUTS_from_m2_f/analysis",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: max(10 * input.size_mb * attempt, 4000),
        time_min = 59
    script:
        "../scripts/STEP7-plot_var_metric_2D.py"
