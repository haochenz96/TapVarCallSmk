rule step3_bcftools_filter_and_merge_sc_m2_vcfs:
    # scatter by sample
    # for each single cell VCF, split multiallelic sites, hard filter (v3) and index; 
    input:
        sc_m2_vcfs_filter_added = get_step2_m2_vcfs_by_sample,
    output:
        merged_prev_filtered_vcf = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
        merged_prev_filtered_vcf_stats = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz.stats",
    params:
        # @HZ 04/01/2022
        # need to add quotation marks around these filters in the shell command, because snakemake will input these params.X as string without quotation marks
        sc_vcf_filters = parse_bcf_filters(config['bcftools']['pass1_sc_vcf_filters']),
        merged_vcf_filters = parse_bcf_filters(config['bcftools']['pass1_merged_vcf_filters']), 
    conda:
        "../envs/bcftools.yaml"
    retries: 1
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    log: 
        "{sample_name}/logs/{sample_name}.STEP3.merge_sc_m2_calls.log"
    shell:
        """
        echo "----------" >> {log}
        echo "[`date`]" >> {log}
        echo "[INFO] [STEP2] finished sc_mutect2_call for sample {wildcards.sample_name}" >> {log}

        # @part1
        mkdir -p {wildcards.sample_name}/3-bcf_filter/filtered_vcfs
        for sc_m2_vcf in {input.sc_m2_vcfs_filter_added}; do 
            
            out_name=$(basename ${{sc_m2_vcf}})

            bcftools norm -m- ${{sc_m2_vcf}} | \
            bcftools filter -i '{params.sc_vcf_filters}' | \
            bgzip -c > \
            {wildcards.sample_name}/3-bcf_filter/filtered_vcfs/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

            tabix {wildcards.sample_name}/3-bcf_filter/filtered_vcfs/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}
            # echo cell -- ${{out_name}} -- norm, filter, zip, index done
        done

        # @part2
        bcftools merge --merge none $(find {wildcards.sample_name}/3-bcf_filter/filtered_vcfs -name '*hard_filtered.vcf.gz') | \
        bcftools sort | \
        bcftools filter -i '{params.merged_vcf_filters}' | \
        bgzip -c > \
        {output.merged_prev_filtered_vcf}

        tabix {output.merged_prev_filtered_vcf}

        bcftools stats {output.merged_prev_filtered_vcf} > {output.merged_prev_filtered_vcf_stats} && \
         echo "[INFO] [STEP3] finished filtering and merging single-cell mutect2(pass1) VCFs sample {wildcards.sample_name}" >> {log}
        """

rule step3_fetch_snv_base_qual_infos:
    input:
        sc_m2_vcfs_filter_added = get_step2_m2_vcfs_by_sample,
    output:
        temp_output_dir = temp(directory("{sample_name}/3-bcf_filter/__bq")),
        merged_snv_bq_vcf = "{sample_name}/3-bcf_filter/merged_bq_info/{sample_name}-bq_merged.sc_prev.csv",
    params:
        sample_name = "{sample_name}",
        output_dir = "{sample_name}/3-bcf_filter/merged_bq_info",
        # sc_filter = 'FILTER~"base_qual"',
    conda:
        "../envs/bcftools.yaml",
    retries: 1
    threads: lambda wildcards, attempt: attempt * 8,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    log: 
        "{sample_name}/logs/{sample_name}.STEP3.merge_sc_m2_calls.log"
    shell:
        """
        echo "+++++++++++++++++++++++++++++++" >> {log}
        echo "----- [time stamp] ----- `date`" >> {log}
        echo "started job: <STEP3> fetch SNV base-quality sc-prev" >> {log}


        # ----- normalize, hard filtering
        # temp_output_dir={output.temp_output_dir}
        mkdir -p {output.temp_output_dir}/__filtered_sc_vcfs

        for sc_m2_vcf in {input.sc_m2_vcfs_filter_added}
        do 
            out_name=$(basename $sc_m2_vcf)
            # echo "processing $out_name"
            
            normed_filtered_sc_vcf={output.temp_output_dir}/__filtered_sc_vcfs/${{out_name/somatic_m2_filter_added.vcf.gz/.normed.bq.vcf.gz}}

            # split multiallelic sites into biallelic records (-m-)
            # ***** CRITICALLY *****: only include SNVs that have BAD base_qual
            bcftools norm -m- $sc_m2_vcf | \
            bcftools filter \
                -i 'FILTER~"base_qual"' \
                -Oz \
                -o $normed_filtered_sc_vcf && \
            tabix $normed_filtered_sc_vcf 
            echo cell -- $out_name -- norm, filter, zip, index done >> {log}
        done

        if [ $? -eq 0 ]; then echo '[INFO] finished normalizing BQ-filtered VCFs' >> {log}; fi

        bcftools merge $(find {output.temp_output_dir}/__filtered_sc_vcfs -name '*.normed.bq.vcf.gz')\
            --merge none | \
        bcftools sort | \
        bgzip -c \
            > {params.output_dir}/{params.sample_name}-bq_merged.vcf.gz && \
        tabix {params.output_dir}/{params.sample_name}-bq_merged.vcf.gz && \
        echo '[INFO] finished merging VCFs' >> {log}

        # ----- get vcf stats
        bcftools stats {params.output_dir}/{params.sample_name}-bq_merged.vcf.gz > {params.output_dir}/{params.sample_name}-bq_merged.vcf.gz.stats
        # get sc-prev for each SNV
        bcftools query {params.output_dir}/{params.sample_name}-bq_merged.vcf.gz \
            -f '%CHROM:%POS:%REF/%ALT,%N_PASS(GT!="mis")\n' \
            > {params.output_dir}/{params.sample_name}-bq_merged.sc_prev.csv
        """