###########
# ### STEP1 get step1 sc_bams output
###########
def get_step1_sc_bams(sample_names, sample_barcode_maps):
    # for the main Snakefile as target output
    # gets all sc_bams for ALL sample
    out = []
    for sample_i in sample_names:
        for cell_num_index in list(sample_barcode_maps[sample_i].keys()):
            out.append(f"{sample_i}/1-sc_bams/{sample_i}_{cell_num_index}.bam")
    return out

def get_step1_sc_bams_by_sample(wildcards):
    # for Snakemake rule use
    # gets all sc_bams for A GIVEN sample
    sample_i = wildcards.sample_name
    out = []
    for cell_num_index in list(sample_barcode_maps[sample_i].keys()):
        out.append(f"{sample_i}/1-sc_bams/{sample_i}_{cell_num_index}.bam")
    return out

###########
# ### STEP2 get the single-cell "filter_added.vcf.gz" outputs
###########
def get_step2_filter_added_mutect2_vcfs(sample_names, sample_barcode_maps):
    # for the main Snakefile as TARGET OUTPUT
    out = []
    for sample_i in sample_names:
        for cell_num_index in sample_barcode_maps[sample_i].keys():
            out.append(f"{sample_i}/2-sc_mutect2_call/m2_sc_vcfs_filter_added/{sample_i}_{cell_num_index}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step2_m2_vcfs_by_sample(wildcards):
    # for Snakemake rule use only
    # gets all step2 mutect2 output vcf for A GIVEN sample
    # as input for step2 bcftools filter
    out = []
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/2-sc_mutect2_call/m2_sc_vcfs_filter_added/{wildcards.sample_name}_{cell_num_index}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step2_strelka2_vcfs(sample_names, sample_barcode_maps):
    # for the main Snakefile as target output
    out = []
    for sample_i in sample_names:
        for cell_num_index in sample_barcode_maps[sample_i].keys():
            #"{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/results/variants/variants.vcf.gz",
            out.append(f"{sample_i}/sc_strelka2_call/{sample_i}_{cell_num_index}/results/variants/variants.vcf.gz")
    return out

###########
# ### STEP3 ----- single-cell Mutect2 filtering; merge single cell VCFs
###########
def get_step2_m2_filtered_vcfs_by_sample(wildcards):
    # for Snakemake rule use only
    # gets all step2 bcftools-filtered vcf for A GIVEN sample
    # as input for step3 bcftools merge
    out = []
    for cell_barcode in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/3-bcf_filter/filtered_vcfs/{wildcards.sample_name}_{cell_barcode}_hard_filtered.vcf.gz")
    return out


###########
# ### STEP4 ----- single-cell genotying of SNVs in q_vcf from above
###########
def get_step4_m2_f_vcfs(wildcards):
    out = []
    for cell_barcode in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/4-sc_mutect_f_call/m2_sc_vcfs/{wildcards.sample_name}_{cell_barcode}_somatic_m2.vcf.gz")
    return out

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