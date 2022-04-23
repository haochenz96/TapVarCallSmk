# @HZ 03/09/2022
# prepare all inputs for running custom single-cell variant calling pipeline
# - 

import pandas as pd
from snakemake import load_configfile
from pathlib import Path
import os
import pysam
import glob

config = load_configfile('configs/custom_tap_pipeline_hg19-b37.yaml')

# ----- get sample and reference info -----
cohort_name = config['sample_info']['cohort_name']
sample_names = config['sample_info']['sample_names']
genome_version = config['reference_info']['genome_version']

# ----- get directory info ----------------
top_dir = Path(config['top_dir'])
if config['working_dir'] == "default":
    print(f"[INFO] working_dir is set to default: {config['top_dir']}/{cohort_name}")
    working_dir = top_dir / cohort_name
else:
    print(f"[INFO] working_dir specified to be: {config['working_dir']}")
    working_dir = topdir / config['working_dir']

scripts_dir = Path(config['scripts_dir'])

#print(f'[INFO] wd --- {working_dir}')
# example: /home/zhangh5/work/Tapestri_batch1/pipeline_results_Mutect2/Caitlin/M04-3_hg19-b37

###########
# ### STEP1 ----- Tapestri-pipeline
###########

# part_1_output = working_dir / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_names}_{genome_version}.tube1.cells.bam'
#print(f'[INFO] part-1 output --- {part_1_output}')

###########
# ### STEP2 ----- single-cell Mutect2 calling
###########

# beginning input: combined_cells.bam
# combined_cells_bam = config['combined_cells.bam']
# if not part_1_output.is_file():
#     print("[WARNING] part 1 output -- combined_cells.bam -- does not exist!")

# get single-cell barcodes for each sample
sample_barcode_maps = {} # <--------------------------------------------- global
for sample_i in sample_names:
    bars_map = {} # <------------------------------------------------------- global
    part_1_output = working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    with pysam.AlignmentFile(part_1_output, "rb") as cells_bam:
        # @HZ 04/19/2022 create numerical index for each barcode
        cell_count = 0
        for i in cells_bam.header['RG']:
            cell_count += 1
            bars_map[f'cell_{cell_count}'] = i['SM'] # <-------------------- numerical index to cell barcode map
            #bars.append(i['SM'])

    sample_barcode_maps[sample_i] = bars_map
print('[INFO] barcode reading done!')

# (1) get step1 sc_bams output
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

# (2) get the single-cell "filter_added.vcf.gz" outputs
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
    for cell_barcode in sample_barcode_mapss[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/3-bcf_filter/filtered_vcfs/{wildcards.sample_name}_{cell_barcode}_hard_filtered.vcf.gz")
    return out


###########
# ### STEP4 ----- single-cell Mutect2 force-calling using q_vcf from above
###########

# def get_step4_outputs(sample_names, sample_barcode_maps):
#     # for the main Snakefile as target output
#     out = []
#     for sample_i in sample_names:
#         for cell_barcode in sample_barcode_maps[sample_i]:
#             out.append(f"{sample_i}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_i}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
#     return out

# def get_step4_m2_f_vcfs(wildcards):
#     out = []
#     for cell_barcode in sample_barcode_maps[wildcards.sample_name]:
#         out.append(f"{wildcards.sample_name}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{wildcards.sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
#     return out

def get_step4_mpileup_vcfs(wildcards):
    out = []
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_num_index}_mpileup.vcf.gz")
    return out