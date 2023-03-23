rule step4_generate_candidate_allele:
    # generate a candidate allele from Q_VCF
    # @HZ 08/28/2022: need to merge multiallelic sites for the next step (mpileup+call) to work
    input:
        Q_VCF = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output: 
    	CANDIDATE_ALLELE = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.atomized.multiallelic.tsv.gz",
        CANDIDATE_ALLELE_for_py = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.atomized.multiallelic.for_py.csv",
    threads: lambda wildcards, attempt: attempt * 2
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/bcftools.yaml"
    params:
        Q_VCF_atomized = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.atomized.vcf.gz",
        Q_VCF_multiallelic = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.atomized.multiallelic.vcf.gz",
        REF_GENOME = config['reference_info']['reference_genome'],
        PANEL_BED = config['reference_info']['panel_insert_file'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    log: 
        "{sample_name}/logs/STEP4-mpileup_genotype.log"
    shell:
        """
        TIMESTAMP=[`date`]
        echo $TIMESTAMP >> {log} 
        echo '--- STEP4 - generating candidate alleles from combined, filtered sc-mutect2 VCFs.' >> {log} 
        
        # @HZ 03/21/2023 atomize MNPs
        {params.BCFTOOLS_EXEC} norm \
            -Oz \
            -o {params.Q_VCF_atomized} \
            --atomize \
            {input.Q_VCF}

        # @HZ 09/18/2022 get a condensed form for Python use later
        {params.BCFTOOLS_EXEC} query \
            -f '%CHROM:%POS:%REF/%ALT\n' \
            {params.Q_VCF_atomized} > \
            {output.CANDIDATE_ALLELE_for_py}

        # @HZ 08/28/2022: merge multiallelic sites 
        {params.BCFTOOLS_EXEC} norm \
            -Oz \
            -o {params.Q_VCF_multiallelic} \
            --multiallelics + \
            {params.Q_VCF_atomized}

        {params.BCFTOOLS_EXEC} query \
            -f'%CHROM\t%POS\t%REF,%ALT\n' \
            {params.Q_VCF_multiallelic} | \
        bgzip -c > {output.CANDIDATE_ALLELE} && \
        tabix -s1 -b2 -e2 {output.CANDIDATE_ALLELE} && \
        echo '--- [STEP4] - finished writing candidate alleles file.' >> {log}

        echo '--- [STEP4] - next, bcftools mpileup with:' >> {log} && \
        echo '--- --- BED FILE: {params.PANEL_BED}' >> {log} && \
        echo '--- --- GENOME FILE: {params.REF_GENOME}' >> {log} 
        """

rule step4_sc_mpileup:
    # scattered by single cell
    # mpileup and call variants 
    input:
        SC_BAM = "{sample_name}/1-sc_bams/{sample_name}_{cell_num_index}.bam",
    	CANDIDATE_ALLELE = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.atomized.multiallelic.tsv.gz",
        # HEADER_FILE = "{sample_name}/4-bcf_genotyping/candidate_alleles.header.txt",
    output:
        # SC_RAW_COUNTS_TXT = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts.txt.gz",
        SC_MPILEUP_VCF = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.vcf.gz",
        SC_MPILEUP_VCF_TBI = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.vcf.gz.tbi",
        SC_MPILEUP_AD_LAYER = "{sample_name}/4-bcf_genotyping/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.AD.csv",
        SC_MPILEUP_DP_LAYER = "{sample_name}/4-bcf_genotyping/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
    params:
        REF_GENOME = config['reference_info']['reference_genome'],
        # PANEL_BED = config['reference_info']['panel_insert_file'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    threads: 2
    resources: 
        mem_mb = lambda wildcards, attempt: min(8000, attempt * 4000),
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # --- 1. mpileup ---
        # @HZ 09/19/2022: need to add `-f [REF_GENOME]` to the norm command for indels to be normalized
        # otherwise the merging step won't consider them as they are not present in the alleles file
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
            -f '%CHROM:%POS,%DP\n' \
            > {output.SC_MPILEUP_DP_LAYER}
        """

# # life-saving rule to be enabled when necessary
# rule extract_sc_mpileup_raw:
#     input:
#         SC_MPILEUP_VCF = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.vcf.gz",
#     output:
#         SC_MPILEUP_DP_LAYER = "{sample_name}/4-bcf_genotyping/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
#     params:
#         BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
#     shell:
#         """
#         {params.BCFTOOLS_EXEC} query \
#             {input.SC_MPILEUP_VCF} \
#             -f '%CHROM:%POS:%REF,%DP\n' \
#             > {output.SC_MPILEUP_DP_LAYER}
#         """

def get_step4_sc_mpileup_vcfs(wildcards):
    out = []
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_num_index}_raw_counts_added.vcf.gz")
    return out

def get_step4_sc_mpileup_AD(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_raw_data/{wildcards.sample_name}_{cell_num_index}_mpileup.normed.AD.csv")
    return out

def get_step4_sc_mpileup_DP(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_raw_data/{wildcards.sample_name}_{cell_num_index}_mpileup.normed.DP.csv")
    return out

rule step4_merge_raw_layers:
    # scatter by sample
    # merge raw counts layers
    input:
        SC_MPILEUP_AD_LAYER = get_step4_sc_mpileup_AD,
        SC_MPILEUP_DP_LAYER = get_step4_sc_mpileup_DP,
        CANDIDATE_ALLELE_for_py = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.atomized.multiallelic.for_py.csv",
    output:
        SC_MPILEUP_AD_LAYER_MERGED = "{sample_name}/4-bcf_genotyping/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "{sample_name}/4-bcf_genotyping/{sample_name}.mpileup.DP.merged.csv",
    params:
        input_dir = "{sample_name}/4-bcf_genotyping/sc_mpileup_raw_data",
        output_dir = "{sample_name}/4-bcf_genotyping",
        MERGE_SCRIPT = config['python_scripts']['MERGE_MPILEUP_RAW_DATA_SCRIPT'],
    log:
        std = "{sample_name}/logs/{sample_name}.STEP4-mpileup_genotype.log",
        err = "{sample_name}/logs/{sample_name}.STEP4-mpileup_genotype.err",
    threads: lambda wildcards, attempt: attempt * 4,
    retries: 2,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = lambda wildcards, attempt: attempt * 59,
    shell:
        """
        python {params.MERGE_SCRIPT} \
            --sample_name {wildcards.sample_name} \
            --input_dir {params.input_dir} \
            --candidate_alleles_df {input.CANDIDATE_ALLELE_for_py} \
            --output_dir {params.output_dir} \
            1> {log.std} 2> {log.err}
        """