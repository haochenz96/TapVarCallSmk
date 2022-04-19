rule bcftools_pass_2_filter_merge_sc_m2_calls:
    # scatter by sample
    # (1) for each single cell VCF, split multiallelic sites, apply only DP filter and index; 
    # (2) merge filtered single cell VCFs; no filter based on mutational prevalence across all cells; create statistics.
    input:
        sc_f_vcfs_filter_added = get_step4_m2_f_vcfs,
    output:
        merged_f_vcf = "{sample_name}/bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        merged_f_vcf_stats = "{sample_name}/bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz.stats",
    params:
        sc_vcf_filters = parse_bcf_filters(config['bcftools']['pass2_sc_vcf_filters']),
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    shell:
        """
        mkdir -p {wildcards.sample_name}/bcf_pass2/filtered_vcfs
        for m2_vcf in {input.sc_f_vcfs_filter_added}; do 
            
            out_name=$(basename ${{m2_vcf}})
            
            
            bcftools norm -m- ${{m2_vcf}} | \
            bcftools filter -i '{params.sc_vcf_filters}' | \
            bgzip -c > \
            {wildcards.sample_name}/bcf_pass2/filtered_vcfs/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}

            tabix {wildcards.sample_name}/bcf_pass2/filtered_vcfs/${{out_name/somatic_m2_filter_added.vcf.gz/hard_filtered.vcf.gz}}
            echo cell -- ${{out_name}} -- norm, filter, zip, index done
        done

        bcftools merge --merge none $(find {wildcards.sample_name}/bcf_pass2/filtered_vcfs -name '*hard_filtered.vcf.gz') | \
        bcftools sort | \
        bgzip -c > \
        {output.merged_f_vcf}

        tabix {output.merged_f_vcf}

        bcftools stats {output.merged_f_vcf} > {output.merged_f_vcf_stats}
        """

rule bcftools_pass2_intersect_f_q_vcfs:
    # intersect q_vcf and f_vcf
    # all samples together
    input:
        f_vcf = "{sample_name}/bcf_pass2/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        q_vcf = "{sample_name}/bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output:
        output_vcf = "{sample_name}/bcf_pass2/{sample_name}-f_q_intersected.vcf.gz",
        output_vcf_stats = "{sample_name}/bcf_pass2/{sample_name}-f_q_intersected.vcf.gz.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        bcftools isec \
            -n=2 \
            -w2 \
            -Oz \
            -o {output.output_vcf} \
            {input.q_vcf} \
            {input.f_vcf}
        
        tabix {output.output_vcf}
        bcftools stats {output.output_vcf} > {output.output_vcf_stats}
        """