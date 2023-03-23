# @HZ 06/20/2022
# script to add CRAVAT annotations to H5 
# inputs:
#   - an H5 file annotated with bulk information
#   - a file with the variants and their CRAVAT annotations
# outputs:
#   - the annotated H5 file

from tea.parse import *
from tea.format import isNaN
from tea.utils import get_simple_timestamp
import mosaic.io as mio
import sys

# # ----- io -----
# SAMPLE_NAME = snakemake.wildcards.sample_name
# # --- inputs
# H5 = snakemake.input.ponv1_filtered_h5
# CRAVAT_OUTPUT = Path(snakemake.input.cravat_output)
# # --- outputs
# ANNOTATED_H5 = snakemake.output.annotated_h5
# CRAVAT_OUTPUT_CLEANED = Path(snakemake.output.cravat_output_cleaned)
# # --- params
# LOG_FILE = Path(snakemake.params.log_file)

# # for testing
# SAMPLE_NAME = 'M13-1'
# H5 = '/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/M13-1_combined_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5'
# __output_dir = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/test')
# __output_dir.mkdir(exist_ok=True, parents=True)
# CRAVAT_OUTPUT = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/annotations/M13-1_combined_CRAVAT_output.txt')
# ANNOTATED_H5 = __output_dir / 'M13-1_combined_annotated.h5'
# CRAVAT_OUTPUT_CLEANED = __output_dir / 'M13-1_combined_CRAVAT_output_cleaned.txt'
# LOG_FILE = __output_dir / 'STEP6_CRAVAT_log.txt'

simple_timestamp = get_simple_timestamp()

with open (LOG_FILE, 'a') as f:
    sys.stdout = f # redirect stdout to log file
    print(f'--- {simple_timestamp} ---')
    # 0. read in H5
    sample_obj = mio.load(H5, name = SAMPLE_NAME, raw = False)

    # 1. read in the CRAVAT report
    with read_until_match(
        in_file = CRAVAT_OUTPUT, 
        match_pattern = (lambda line: "#CRAVAT Report" in line), 
        num_occurences = 2
        ) as f:
        
        f.seek(0) # !!!reset cursor!!!
        cravat_df = pd.read_csv(f, sep='\t', skiprows=4, header=[0,1], low_memory=False)

    # 2. clean and format the CRAVAT report
    cravat_df = clean_and_format_cravat_df(cravat_df)


    # 3. add key information from H5 variant annotations
    for ds in sample_obj.dna.col_attrs.keys():
        if 'bulk' in ds:
            print(f'[INFO] adding bulk comparisons for {ds}')
            if type(sample_obj.dna.col_attrs[ds]) == str:
                sample_obj.dna.col_attrs[ds] = bool(sample_obj.dna.col_attrs[ds])
            cravat_df.loc[:, ('bulk_comparison', ds)] = sample_obj.dna.col_attrs[ds]
        elif 'pon' in ds:
            print(f'[INFO] adding PoN comparisons for {ds}')
            if type(sample_obj.dna.col_attrs[ds]) == str:
                sample_obj.dna.col_attrs[ds] = bool(sample_obj.dna.col_attrs[ds])
            cravat_df.loc[:, ('pon_comparison', ds)] = sample_obj.dna.col_attrs[ds]
        elif ds == 'sc_max_TLOD':
            print(f'[INFO] adding SC max TLOD for {ds}')
            cravat_df.loc[:, ('Tapestri_result', ds)] = sample_obj.dna.col_attrs[ds]
        else:
            continue
    # adjust the order of bulk & pon comparison columns
    for lv1_col_name in ['bulk_comparison', 'pon_comparison']:
        cols = cravat_df.pop(lv1_col_name)
        for i in range(cols.shape[1]):
            cravat_df.insert(i, (lv1_col_name, cols.columns[i]), cols.iloc[:,i])
    print(f'[INFO] pushed forward the columns for {ds}')

    # insert TLOD max value into the DataFrame
    i = np.where(cravat_df.columns == ('Tapestri_result', 'sc_mut_prev'))[0][0] + 1
    tlod_col = cravat_df.pop(('Tapestri_result', 'sc_max_TLOD'))
    cravat_df.insert(int(i), ('Tapestri_result', 'sc_max_TLOD'), tlod_col)
    
    # important: fetch variant annotation before any sorting
    gene_HGVSp = cravat_df.index.map(
        lambda x: 
        cravat_df.loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation', 'Protein Change')] if not cravat_df.loc[x, ('Variant Annotation','Protein Change')]
        else cravat_df.loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation','Sequence Ontology')]
    )
    # from IPython import embed
    # embed()    
    sample_obj.dna.add_col_attr('gene_HGVSp', np.array(gene_HGVSp))

    # sort by Tapestri sc mutational prevalence
    cravat_df.sort_values(by = ('Tapestri_result', 'sc_mut_prev'), ascending=False, inplace=True)
    
    # write to outputs:
    mio.save(sample_obj, ANNOTATED_H5)
    print(f'[INFO] saved annotated H5 file to {ANNOTATED_H5}')
    cravat_df.to_csv(ANNOTATED_H5.with_suffix('.cravat.tsv'), sep='\t', index=False)
    print(f'[INFO] saved cleaned CRAVAT report to {CRAVAT_OUTPUT_CLEANED}')

