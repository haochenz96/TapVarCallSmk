rule patient_step3_write_h5_from_mpileup_raw_data_and_rc:
    # write SNV and CNV matrices.
    input:
        SC_MPILEUP_AD_LAYER_MERGED = "fillout/{sample_name}/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "fillout/{sample_name}/{sample_name}.mpileup.DP.merged.csv",
    output:
        output_h5 = "fillout/{sample_name}/{sample_name}.mpileup.h5",
    params:
        write_h5_script = os.path.join(scripts_dir, "single_sample/STEP5-write_h5_from_raw_AD_DP_and_rc.py")
    log: 
        "fillout/logs/{sample_name}-genotype.log",
    conda:
        "../../envs/mosaic-custom.yaml"
    threads: 4
    # retries: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = 59
    shell:
        """
        python {params.write_h5_script} \
            --sample_name {wildcards.sample_name} \
            --AD_merged_df {input.SC_MPILEUP_AD_LAYER_MERGED} \
            --DP_merged_df {input.SC_MPILEUP_DP_LAYER_MERGED} \
            --output_h5 {output.output_h5} \
            > {log} 2>&1
        """

rule patient_step3_combine_sample_h5s:
    input:
        sample_h5s = expand("fillout/{sample_name}/{sample_name}.mpileup.h5", sample_name = sample_names),
    output:
        merged_h5 = "fillout/{patient_name}.patient_wide.genotyped.h5", 
    params:
        sample_names = sample_names
    log: 
        "fillout/logs/{patient_name}-genotype.log",
    conda:
        "../../envs/mosaic-custom.yaml"
    threads: 4
    # retries: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = 59
    shell:
        """
        echo "----------" >> {log}
        echo "[`date`]" >> {log}
        echo "[INFO] [PATIENT-STEP3] merging genotyped H5s for samples {params.sample_names}" >> {log}
        tapestri h5 merge samples \
            {input.sample_h5s} \
            {output.merged_h5} \
            >> {log} 2>&1
        """