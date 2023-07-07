rule patient_sc_mpileup:
    # scattered by each sample, each single cell
    # mpileup and call variants 
    input:
        SC_BAM = get_individual_sc_bam, # retrieve EACH sc_bam given the wildcards: sample_name and cell_num_index
        CANDIDATE_ALLELE = expand("fillout/input/{patient_name}.snv_union.for_genotyping.vcf.gz", patient_name = patient_name),
    output:
        # SC_RAW_COUNTS_TXT = "fillout/{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts.txt.gz",
        SC_MPILEUP_VCF = "fillout/{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz",
        SC_MPILEUP_VCF_TBI = "fillout/{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz.tbi",
        # SC_MPILEUP_VCF_raw_counts = "fillout/{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts_added.vcf.gz",
        SC_MPILEUP_AD_LAYER = "fillout/{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.AD.csv",
        SC_MPILEUP_DP_LAYER = "fillout/{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
    params:
        REF_GENOME = config['reference_info']['reference_genome'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    threads: 2
    group: "sc_mpileup"
    resources: 
        mem_mb = lambda wildcards, attempt: min(8000, attempt * 4000),
        time_min = lambda wildcards, attempt: attempt * 119,
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        # --- 1. mpileup ---
        # @HZ 09/19/2022: need to add `-f [REF_GENOME]` to the norm command for indels to be normalized
        # otherwise the merging step won't consider them as they are not present in the alleles file

        # @HZ 03/08/2023: added sort step because the candidate alleles file could be unsorted
        {params.BCFTOOLS_EXEC} mpileup \
            {input.SC_BAM} \
            -R {input.CANDIDATE_ALLELE} \
            -f {params.REF_GENOME} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            --max-depth 100000 \
            --max-idepth 100000 \
            -Ou | \
        {params.BCFTOOLS_EXEC} norm \
        	-m- \
            -f {params.REF_GENOME} \
        	-Ou | \
        {params.BCFTOOLS_EXEC} sort \
            -Oz \
            > {output.SC_MPILEUP_VCF} && \
        tabix {output.SC_MPILEUP_VCF}

        # --- 2. extract AD, DP ---
        {params.BCFTOOLS_EXEC} query \
        	{output.SC_MPILEUP_VCF} \
        	-f '%CHROM:%POS:%REF/%ALT,%AD{{1}}\n' \
        	-i'ALT!="<*>"' \
        	> {output.SC_MPILEUP_AD_LAYER}
        {params.BCFTOOLS_EXEC} query \
            {output.SC_MPILEUP_VCF} \
            -f '%CHROM:%POS:%REF,%DP\n' \
            > {output.SC_MPILEUP_DP_LAYER}
        """

# rule extract_sc_mpileup_raw:
#     input:
#         SC_MPILEUP_VCF = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz",
#     output:
#         SC_MPILEUP_DP_LAYER = "{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
#     params:
#         BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
#     shell:
#         """
#         {params.BCFTOOLS_EXEC} query \
#             {input.SC_MPILEUP_VCF} \
#             -f '%CHROM:%POS:%REF,%DP\n' \
#             > {output.SC_MPILEUP_DP_LAYER}
#         """

rule merge_raw_layers:
    # scatter by sample
    # merge raw counts layers
    input:
        SC_MPILEUP_AD_LAYER = get_sc_mpileup_AD,
        SC_MPILEUP_DP_LAYER = get_sc_mpileup_DP,
        CANDIDATE_ALLELE_for_py = expand("fillout/input/{patient_name}.snv_union.multiallelic.for_py.csv", patient_name = patient_name),
    output:
        SC_MPILEUP_AD_LAYER_MERGED = "fillout/{sample_name}/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "fillout/{sample_name}/{sample_name}.mpileup.DP.merged.csv",
    params:
        input_dir = "fillout/{sample_name}/sc_mpileup_raw_data",
        output_dir = "fillout/{sample_name}",
        MERGE_SCRIPT = os.path.join(scripts_dir, "single_sample/STEP4-merge_sc_AD_DP.py")
    log:
        "fillout/logs/{sample_name}-genotype.log",
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 59,
    shell:
        """
        python {params.MERGE_SCRIPT} \
            --sample_name {wildcards.sample_name} \
            --input_dir {params.input_dir} \
            --candidate_alleles_df {input.CANDIDATE_ALLELE_for_py} \
            --output_dir {params.output_dir} \
            > {log} 2>&1
        """