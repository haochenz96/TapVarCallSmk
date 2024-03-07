from pathlib import Path
import sys
import glob
import re

# ----- input functions -----
def get_candidate_alleles(patient_name):
    if config['patient_info']['candidate_alleles'] is None:
        return f"fillout/input/{patient_name}.snv_union.for_genotyping.vcf.gz"
    else:
        return config['patient_info']['candidate_alleles']

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
#         patient_wide_workflow_outputs,

