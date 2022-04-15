from missionbio.h5.create import create_cnv_assay, create_dna_assay
from missionbio.h5.data import H5Writer
from missionbio.h5.constants import BARCODE, CHROM, ID, POS
import missionbio.mosaic.io as mio
import missionbio.mosaic.utils
import allel
import pandas as pd
import numpy as np
import click
import json
import os

###############################
# part 0 ----- parse inputs
###############################

# @click.option('--sample_name', required=True, type=str)
# @click.option('--metadata', required=True, type=str)
# @click.option('--all_cell_vcf', required=True, type=str)
# @click.option('--amplicon_file', required=True, type=str)
# @click.option('--read_count_tsv', required=True, type=str)
# @click.option('--output_dir', required=True, type=str)

# metadata = json.loads(metadata)

# get variables from snakemake
sample_name = snakemake.wildcards.sample_name

#metadata = snakemake.params.metadata_json
metadata = json.loads(snakemake.params.metadata_json)
for i in metadata:
    if isinstance(metadata[i], dict):
        metadata[i] = json.dumps(metadata[i])
print(metadata)

all_cell_vcf = snakemake.input.output_vcf
amplicons_file = snakemake.config['reference_info']['panel_amplicon_file']
read_counts_tsv = snakemake.input.read_counts_tsv
output_h5 = snakemake.output.output_h5

print(f'[INFO] ----- output file: {output_h5}')

# output_dir = output_h5.split('/')[0]

# print(f'[INFO] ----- output directory: {output_dir}')
# try:
#     os.mkdir(output_dir)
# except OSError as error:
#     print(f'[WARNING] ----- {error}')    

###############################
# part 1 ----- create DNA assay
###############################

dna = create_dna_assay(all_cell_vcf, metadata)

###############################
# part 2 ----- create CNV assay
###############################

def add_amplicon_metadata(cnv_assay, amplicons):
    ca = cnv_assay.col_attrs
    chrom = np.full((cnv_assay.shape[1],), "", dtype=object)
    start_pos = np.full((cnv_assay.shape[1],), 0, dtype=int)
    end_pos = np.full((cnv_assay.shape[1],), 0, dtype=int)

    for _, amplicon in amplicons.iterrows():
        matching_ids = ca[ID] == amplicon['amplicon']
        chrom[matching_ids] = amplicon['CHROM'].strip("chrCHR")
        start_pos[matching_ids] = amplicon['start_pos']
        end_pos[matching_ids] = amplicon['end_pos']

    cnv_assay.add_col_attr('CHROM', chrom)
    cnv_assay.add_col_attr('start_pos', start_pos)
    cnv_assay.add_col_attr('end_pos', end_pos)

cnv = create_cnv_assay(read_counts_tsv, metadata)

amplicons = pd.read_csv(
            amplicons_file,
            sep="\t",
            lineterminator="\n",
            names=['CHROM', 'start_pos', 'end_pos', 'amplicon'],
            dtype={'start_pos': int, 'end_pos': int},
)
add_amplicon_metadata(cnv, amplicons)

###############################
# part 3 ----- add both to H5
###############################

assays = [dna, cnv]

# !rm $output_h5
with H5Writer(output_h5) as writer:
    for assay in assays:
        writer.write(assay)