from pathlib import Path
import sys
import glob
import re

# ----- input functions -----
def get_individual_sc_bam(wildcards):
    # get one single cell's bam file
    return f"{wildcards.sample_name}/1-sc_bams/{wildcards.sample_name}_{wildcards.cell_num_index}.bam"

def get_sc_mpileup_vcfs(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"fillout/{wildcards.sample_name}/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_num_index}_mpileup.normed.vcf.gz")
    print(f"for {wildcards.sample_name}, {len(out)} single cells' VCFs found!")
    return out

def get_sc_mpileup_AD(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"fillout/{wildcards.sample_name}/sc_mpileup_raw_data/{wildcards.sample_name}_{cell_num_index}_mpileup.normed.AD.csv")
    return out

def get_sc_mpileup_DP(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"fillout/{wildcards.sample_name}/sc_mpileup_raw_data/{wildcards.sample_name}_{cell_num_index}_mpileup.normed.DP.csv")
    return out

# # ===== target ======
# rule patient_all_outputs:
#     input:
#         # "input/candidate_alleles.multiallelic.for_py.csv"
#         expand(
#             "fillout/{sample_name}/{sample_name}.mpileup.DP.merged.csv",
#             sample_name=sample_names
#         ),
#         # merged_mpileup_vcf = expand(
#         #     "fillout/{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz",
#         #     sample_name = sample_names,
#         # ),
#         # merged_mpileup_vcf_AD_py = expand(
#         #     "fillout/{sample_name}/{sample_name}-genotyped_combined_AD_for_py.txt",
#         #     sample_name = sample_names
#         # ),
# # ===================

