# @HZ 03/09/2022
# prepare all inputs for running custom single-cell variant calling pipeline
# - 

from snakemake import load_configfile
from pathlib import Path
import os
import pysam
import glob
config = load_configfile('configs/custom_sc_varcall.yaml')

# ----- get sample and reference info -----
cohort_name = config['sample_info']['cohort_name']
sample_name = config['sample_info']['sample_name']
genome_version = config['reference_info']['genome_version']

# ----- get directory info ----------------
top_dir = Path(config['top_dir'])
working_dir = top_dir / cohort_name / f"{sample_name}_{genome_version}" 
scripts_dir = Path(config['scripts_dir'])

#print(f'[INFO] wd --- {working_dir}')
# example: /home/zhangh5/work/Tapestri_batch1/pipeline_results_Mutect2/Caitlin/M04-3_hg19-b37

###########
# ### PART1 ----- Tapestri-pipeline
###########

part_1_output = working_dir / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_name}_{genome_version}.tube1.cells.bam'
#print(f'[INFO] part-1 output --- {part_1_output}')

###########
# ### PART2 ----- single-cell Mutect2 calling
###########

# beginning input: combined_cells.bam
# combined_cells_bam = config['combined_cells.bam']
if not part_1_output.is_file():
    print("[WARNING] part 1 output -- combined_cells.bam -- does not exist!")

# (1) get single-cell barcodes
bars = [] # <---------------------- global
with pysam.AlignmentFile(part_1_output, "rb") as cells_bam:
    for i in cells_bam.header['RG']:
        bars.append(i['SM'])
#print(bars)

# (2) 
def get_single_cell_name_mapping(sc_bam_dir, out_file):

    
    print(f'[DEBUG] --- {out_file}')
    with open(out_file, 'w') as out:
        
        count = 0
        for sc_bam in glob.glob(f'{sc_bam_dir}/*.bam'):
            count += 1
            out.write(sc_bam + '\t' + f'cell_{count}')
            out.write('\n')

###########
# ### PART3 ----- single-cell Mutect2 filtering; merge single cell VCFs
###########

def parse_bcf_filters(bcf_filters):
    '''
    input:
        - bcf_filters: list of bcftools filters to use
    output: 
        - string: e.g. 'FILTER!~"base_qual" & FILTER!~"low_allele_frac" & FILTER!~"weak_evidence" & FILTER!~"slippage" & FILTER!~"multiallelic" & FILTER!~"clustered_events" & FMT/DP>4 & AF>0.2'
    '''
    