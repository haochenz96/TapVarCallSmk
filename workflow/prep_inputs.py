# @HZ 03/09/2022
# prepare all inputs for running custom single-cell variant calling pipeline
# - 

from snakemake import load_configfile
from pathlib import Path
import os
import pysam
import glob

config = load_configfile('configs/custom_tap_pipeline_hg19-ucsc.yaml')

# ----- get sample and reference info -----
cohort_name = config['sample_info']['cohort_name']
sample_names = config['sample_info']['sample_names']
genome_version = config['reference_info']['genome_version']

# ----- get directory info ----------------
top_dir = Path(config['top_dir'])
working_dir = top_dir / cohort_name
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

# (1) get single-cell barcodes for each sample
sample_barcode_map = {} # <--------------------------------------------- global
for sample_i in sample_names:
    bars = [] # <------------------------------------------------------- global
    part_1_output = working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    with pysam.AlignmentFile(part_1_output, "rb") as cells_bam:
        for i in cells_bam.header['RG']:
            bars.append(i['SM'])
    sample_barcode_map[sample_i] = bars
print('[INFO] ----- barcode reading done')

# (2) get the single-cell "filter_added.vcf.gz" outputs
def get_step2_outputs(sample_names, sample_barcode_map):
    # for the main Snakefile as target output
    out = []
    for sample_i in sample_names:
        for cell_barcode in sample_barcode_map[sample_i]:
            out.append(f"{sample_i}/mutect2_sc_pass1/m2_sc_vcfs_filter_added/{sample_i}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step2_m2_vcfs(wildcards):
    # for Snakemake rule use only
    # as input for snakemake's step3
    out = []
    for cell_barcode in sample_barcode_map[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/mutect2_sc_pass1/m2_sc_vcfs_filter_added/{wildcards.sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
    return out



###########
# ### STEP3 ----- single-cell Mutect2 filtering; merge single cell VCFs
###########



###########
# ### STEP4 ----- single-cell Mutect2 force-calling using q_vcf from above
###########

def get_step4_outputs(sample_names, sample_barcode_map):
    out = []
    for sample_i in sample_names:
        for cell_barcode in sample_barcode_map[sample_i]:
            out.append(f"{sample_i}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_i}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step4_m2_f_vcfs(wildcards):
    out = []
    for cell_barcode in sample_barcode_map[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{wildcards.sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step4_mpileup_vcfs(wildcards):
    out = []
    for cell_barcode in sample_barcode_map[wildcards.sample_name]:

        #"{sample_name}/bcf_mpileup/sc_mpileup_vcfs/{sample_name}_{cell_barcode}_mpileup.vcf.gz",
        out.append(f"{wildcards.sample_name}/bcf_mpileup/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_barcode}_mpileup.vcf.gz")
    return out