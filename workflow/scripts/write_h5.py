from missionbio.h5.create import create_cnv_assay, create_dna_assay
from missionbio.h5.data import H5Writer

import missionbio.mosaic.io as mio
import missionbio.mosaic.utils
import pandas as pd

###############################
# part 1 ----- create DNA assay
###############################

vcf_file = 

metadata = {
    'sample_name': sample_name,
    'panel_name': 'iacobuzc1440'
}

dna = create_dna_assay(str(vcf_file), metadata)

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

# read in read count and create CNV assay
amplicons_file = '/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch1/panel/1440/1440.amplicons' # <----
read_counts_file = data_dir / 'RSX381_hg19-b37.tube1.barcode.cell.distribution.tsv' # <----

cnv = create_cnv_assay(read_counts_file, metadata)

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
output_h5 = wd / f'{sample_name}_m2_f.dna.h5'

!rm $output_h5
with H5Writer(output_h5) as writer:
    for assay in assays:
        writer.write(assay)