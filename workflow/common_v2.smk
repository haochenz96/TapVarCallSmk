# @HZ 03/09/2022
# prepare all inputs for running custom single-cell variant calling pipeline
# @HZ 10/01/2022
# V2: allow for custom barcode file

import json
import logging
import shutil, os
from datetime import datetime
import pytz as tz
from snakemake import load_configfile
from pathlib import Path
import pysam
import pandas as pd
import sys, os

# ----- specify directories -----
patient_name = config['patient_info']['patient_name']
sample_names = config['patient_info']['sample_names']
top_dir = Path(config['top_dir'])
if config['working_dir'] == "default":
    print(f"[INFO] working_dir is set to default: {config['top_dir']}/{patient_name}")
    working_dir = top_dir / patient_name
else:
    print(f"[INFO] working_dir specified to be: {config['working_dir']}")
    working_dir = Path(config['working_dir'])
scripts_dir = config['scripts_dir']
print(f"[INFO] scripts dir: {scripts_dir}")

# ----- get current date and time
eastern = tz.timezone('US/Eastern') # <------ uses US eastern time by default
now = datetime.now(tz.utc).astimezone(eastern) 
timestamp = now.strftime("%a %m %d %H:%M:%S %Z %Y")
date_simple = now.strftime("%Y-%m-%d")

def parse_for_shell_args(input_dict):
    '''
    General function to parse a python dictionary into bash arguments
    '''
    args = ""
    for key, value in input_dict.items():
        args += f"-{key} {value} "

    return args

# ----- get individual samples' inputs (BAMs, barcode maps, read counts) -----
sample_names = config['patient_info']['sample_names']

# check BAM and TSV exist for each sample
assert all([os.path.exists(config['sample_info']['sample_bams'][sample_i]) for sample_i in sample_names]), "Must specify input BAMs for each sample in the config file."
assert all([os.path.exists(config['sample_info']['sample_rc_tsvs'][sample_i]) for sample_i in sample_names]), "Must specify input read count TSVs for each sample in the config file."

for sample_i in sample_names:
    (working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam').mkdir(parents=True, exist_ok=True)
    (working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv').mkdir(parents=True, exist_ok=True)
    if not os.path.lexists(str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam')):
        os.symlink(
            config['sample_info']['sample_bams'][sample_i], 
            str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam')
            )
        print(f"Symlinked {config['sample_info']['sample_bams'][sample_i]} to {str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam')}")
    if not os.path.lexists(str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv' / f'{sample_i}.tube1.barcode.cell.distribution.tsv')):
        os.symlink(
            config['sample_info']['sample_rc_tsvs'][sample_i], 
            str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv' / f'{sample_i}.tube1.barcode.cell.distribution.tsv')
            )
        print(f"Symlinked {config['sample_info']['sample_rc_tsvs'][sample_i]} to {str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv' / f'{sample_i}.tube1.barcode.cell.distribution.tsv')}")

# ----- Fetch all single-cell barcodes from each sample's combined BAM ----- 
def fetch_sample_barcode_map(sample_i, bam_file = None, wd = None):
    '''
    Try to fetch the sample barcode maps for valid sc-BAMs to use
    
    '''
    if bam_file is None:
        if wd is None:
            raise(ValueError('Must provide either a BAM file or a working directory to infer the BAM file!'))
        else:
            bam_file = wd / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    
    # write the barcode map file to this directory
    bars_map_file = wd / sample_i / 'reference' / f'{sample_i}.barcode_map.txt'
    if not bars_map_file.is_file():
        print(f'[WARNING] --- {sample_i} barcode map file not found. Fetching barcode map from BAM file.')
        (wd / sample_i / 'reference').mkdir(exist_ok=True, parents=True)
        bars_map = {}
        with pysam.AlignmentFile(bam_file, "rb") as cells_bam:
            # @HZ 04/19/2022 create numerical index for each barcode
            cell_count = 0
            for i in cells_bam.header['RG']:
                cell_count += 1
                bars_map[f'cell_{cell_count}'] = i['SM'] # <-------------------- numerical index to cell barcode map
                #bars.append(i['SM'])

            bars_map_df = pd.DataFrame.from_dict(bars_map, orient='index', columns = ['cell_barcode'])
            bars_map_df.index.name = 'num_index'
            bars_map_df.to_csv(str(bars_map_file), sep='\t', index=True)
        bars_map = bars_map_df.to_dict()['cell_barcode']
    else:
        bars_map_df = pd.read_csv(bars_map_file, sep='\t', index_col=0)
        # sanity check barcode map:
        with pysam.AlignmentFile(bam_file, "rb") as cells_bam:
            barcodes_from_BAM = pd.DataFrame(cells_bam.header['RG'])['SM']

            if not bars_map_df['cell_barcode'].isin(barcodes_from_BAM).all():
                #raise KeyError(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                sys.exit(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
            elif bars_map_df.shape[0] != barcodes_from_BAM.shape[0]:
                sys.exit(f'[ERROR] --- {sample_i} barcode map file length does not match number of barcodes from BAM! Abort.')
            else:
                print(f'[INFO] --- {sample_i} barcode map file is OK.')
                
    bars_map = bars_map_df.to_dict()['cell_barcode']
    return bars_map

sample_barcode_maps = {} # patient-wide, each sample has a barcode map

if config['sample_info']['sample_barcode_maps'] is not None:
    # sanity check:
    assert sample_names == list(config['sample_info']['sample_barcode_maps'].keys()), f"[WARNING] Sample names in sample_bams -- {sample_names} -- and sample_barcodes -- {list(config['sample_info']['sample_barcode_maps'].keys())} -- do not match!"
    sample_barcode_maps = config['sample_info']['sample_barcode_maps']

    for sample_i, map_i in sample_barcode_maps.items():
        if map_i is None:
            print(f'[WARNING] No barcode map provided for sample {sample_i}! Will use all barcodes in provided BAM file.')
            sample_barcode_maps[sample_i] = fetch_sample_barcode_map(
                sample_i,
                config['sample_info']['sample_bams'][sample_i], 
                working_dir
                )
        else:
            bars_map_file_copy = working_dir / sample_i / 'reference' / f'{sample_i}.barcode_map.txt'
            (working_dir / sample_i / 'reference').mkdir(exist_ok=True, parents=True)
            try:
                shutil.copyfile(map_i, bars_map_file_copy)
            except shutil.SameFileError:
                print(f'[WARNING] --- {sample_i} barcode map file is the same as the one in the working directory!')
            bars_map_df = pd.read_csv(bars_map_file_copy, sep='\t', index_col=0)
            # sanity check barcode map:
            with pysam.AlignmentFile(
                config['sample_info']['sample_bams'][sample_i],
                "rb"
                ) as cells_bam:
                barcodes_from_BAM = pd.DataFrame(cells_bam.header['RG'])['SM']
                if not bars_map_df['cell_barcode'].isin(barcodes_from_BAM).all():
                    #raise KeyError(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                    sys.exit(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                # elif bars_map_df.shape[0] != barcodes_from_BAM.shape[0]:
                #     sys.exit(f'[ERROR] --- {sample_i} barcode map file length does not match number of barcodes from BAM! Abort.')
                else:
                    print(f'[INFO] --- {sample_i} barcode map file is OK. Reading in {bars_map_df.shape[0]} barcodes.')
                    
            bars_map = bars_map_df.to_dict()['cell_barcode']
            sample_barcode_maps[sample_i] = bars_map
else:
    print('[INFO] No barcode map provided for any sample! Will use all barcodes in provided BAM file.')
    for sample_i in sample_names:
        sample_barcode_maps[sample_i] = fetch_sample_barcode_map(sample_i, None, working_dir)

print('[INFO] barcode reading done!')

# # ===== workflow info =====
# # ----- single_sample workflow ------
# if config['single_sample']['run']:
#     include: "single_sample_main.smk"
#     single_sample_steps = config['single_sample']
#     single_sample_outputs = []
#     if single_sample_steps['1_split_sc_bams']:
#         include: "rules/single_sample/1-split_sc_bams.smk"
#     if single_sample_steps['2_sc_mutect2_call']:
#         include: "rules/single_sample/2-sc_mutect2_call.smk"
#     if single_sample_steps['3_filter_merge_sc_m2_call']:
#         include: "rules/single_sample/3-filter_merge_sc_m2_call.smk"
#     if single_sample_steps['4_sc_mpileup']:
#         include: "rules/single_sample/4-sc_mpileup.smk"
#     if single_sample_steps['5_write_h5']:
#         include: "rules/single_sample/5-write_h5.smk"
# else:
#     print("[INFO] Skipping single_sample VarCall workflow...")

# # ----- patient_wide workflow ------
# if config['patient_wide']['run']:
#     include: "patient_wide_main.smk"
#     patient_wide_steps = config['patient_wide']
#     if patient_wide_steps['1_get_candidate_alleles']:
#         include: "rules/patient_wide/1-get_candidate_alleles.smk"
#     if patient_wide_steps['2_sc_mpileup']:
#         include: "rules/patient_wide/2-sc_mpileup.smk"
#     if patient_wide_steps['3_write_h5']:
#         include: "rules/patient_wide/3-write_h5.smk"
# else:
#     print("[INFO] Skipping patient-wide VarCall workflow...")