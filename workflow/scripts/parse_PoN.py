import pandas as pd
import pathlib
import re
from functools import reduce

pon_version = 'v3_3sc_5patients'

bcf_isec_dir = pathlib.Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/PoN_assembly/pop_prev_filtered_individual_vcfs/')

# (1) get sample_mapping
bcf_isec_README = bcf_isec_dir / 'README.txt'
sample_map_df = pd.read_csv(
    bcf_isec_README, 
    sep='\t', 
    skiprows=4, 
    names = ['to_be_renamed', 'mid','sc_prev_filtered']
    )
sample_map_df['sample_name'] = sample_map_df['sc_prev_filtered'].str.split('/').str[-1].str.split('-sc_prev').str[:1].str[0]
sample_map_df['sample_idx'] = sample_map_df['to_be_renamed'].str.extract(r'(\d{4})').astype(int)

# (2) read in each individual filterd normal VCF
pop_prev_file = bcf_isec_dir / 'sites.txt'
pop_prev_df = pd.read_csv(
    pop_prev_file, 
    names = ['CHROM', 'POS', 'REF',	'ALT', 'pop_appearance']
    )
ind_filtered_vcfs = list(bcf_isec_dir.glob('*mut_prev.txt'))
ind_filtered_vcfs_map = {}
for ind_filtered_vcf_i in ind_filtered_vcfs:
    regex = '\d{4}'
    vcf_idx = int(re.findall(regex, str(ind_filtered_vcf_i))[0]) # get the index of the VCF from file name
    sample_name = sample_map_df.loc[sample_map_df['sample_idx'] == vcf_idx, 'sample_name'].values[0]
    ind_filtered_vcf_i_df = pd.read_csv(
        ind_filtered_vcf_i, 
        names = ['CHROM', 'POS', 'REF',	'ALT', f'{sample_name}-sc_prev'],
        skiprows=1,
        sep='\t'
        )
    if not ind_filtered_vcf_i_df['CHROM'].astype(str).str.startswith('chr').any():
        ind_filtered_vcf_i_df['CHROM'] = 'chr' + ind_filtered_vcf_i_df['CHROM'].astype(str) # unify chromosome annotation
    ind_filtered_vcf_i_df['condensed_format'] = ind_filtered_vcf_i_df['CHROM'].astype(str) + ':' \
            + ind_filtered_vcf_i_df['POS'].astype(str) + ':' \
            + ind_filtered_vcf_i_df['REF'].astype(str) + '/' + ind_filtered_vcf_i_df['ALT'].astype(str)
    ind_filtered_vcf_i_df.set_index('condensed_format', inplace=True)
    ind_filtered_vcf_i_df.drop(columns=['CHROM', 'POS', 'REF',	'ALT'], inplace=True)
    # save to map
    ind_filtered_vcfs_map[sample_name] = ind_filtered_vcf_i_df

# (3)
df_merged = reduce(lambda left, right: pd.merge(left, right, on='condensed_format',
                                            how='outer'), ind_filtered_vcfs_map.values()).fillna(0)
df_merged.to_csv(
    f'/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/PoN_assembly/PoN_{pon_version}.txt', 
    sep='\t')

    


