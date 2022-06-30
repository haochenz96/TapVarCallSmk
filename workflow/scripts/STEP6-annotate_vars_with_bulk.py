# @HZ 06/28/2022
# annotate a variant list with matched bulk info
# integrated into main Snakemake workflow


import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import yaml
import mosaic.io as mio
from tea.utils import read_vcf_to_df
from tea.format import *
from IPython import embed
import json
from typing import Tuple
import re
import sys

# ----- fetch inputs -----
H5 = snakemake.input.input_h5
SAMPLE_NAME = snakemake.wildcards.sample_name
BULK_INFO_YAML= snakemake.input.bulk_info_yaml
OUTPUT_DIR = Path(f'{SAMPLE_NAME}/OUTPUTS_from_m2_f')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
LOG_FILE = snakemake.params.log_file

with open(LOG_FILE, "a") as f:
    sys.stdout = f

    # (0) ----- preprocess single-cell SNV matrix
    sample_obj = mio.load(H5, SAMPLE_NAME, raw = False)
    TLOD_DF = sample_obj.dna.get_attribute('TLOD', constraint='row+col').T.astype(float)
    print('-----------------------------------------------------')
    print(f'TLOD_DF: {TLOD_DF.head(5)}')
    print('-----------------------------------------------------')
    VAR_DF = pd.DataFrame(index=TLOD_DF.index)
    VAR_DF['sc_max_TLOD'] = TLOD_DF.max(axis=1, skipna=True)
    VAR_DF['sc_mean_TLOD'] = TLOD_DF.mean(axis=1, skipna=True)



    # (1) ----- read in bulk info -----
    with open(BULK_INFO_YAML, "r") as f:
        bulk_info = yaml.safe_load(f)

    case_name = bulk_info['case_name']
    # sample_name = bulk_info["sample_name"]

    bulk_dfs = {}
    # -- a. matched bulk tumor; Broad WES PoN
    for ds in ['matched_bulk_tumor', 'broad_wes_pon']:
        if bulk_info[ds] is None:
            print(f"[WARNING] --- {ds} is not given. Skipping.")
            bulk_dfs[ds] = None
            continue
        if not isinstance(bulk_info[ds], str):
            print(f"[ERROR] --- {ds} is not given as a single path!")
        if bulk_info[ds].endswith(('txt', 'maf')):
            bulk_dfs[ds] = pd.read_csv(bulk_info[ds], sep="\t")
        elif bulk_info[ds].endswith(('vcf', 'vcf.gz')):
            bulk_dfs[ds] = read_vcf_to_df(bulk_info[ds])
        else:
            print(f'[ERROR] --- {bulk_info[ds]} is not given as a TXT or MAF file!')

        # get condensed variant format
        if bulk_dfs[ds].index.name not in ('Tapestri_format', 'condensed_format'):
            if ('Tapestri_format' or 'condensed_format') not in bulk_dfs[ds].columns:
                # embed()
                raise ValueError(f'[ERROR] --- {ds} does not have Tapestri_format or condensed_format column!')
                # try:
                #     bulk_dfs[ds]['Tapestri_format'] = bulk_dfs[ds]['CHROM'] + ":" + bulk_dfs[ds]['POS'] + ":" + bulk_dfs[ds]['REF'] + "/" + bulk_dfs[ds]['ALT']
                    
                # except KeyError:
                #     embed()
                #     print(f'[ERROR] --- the columns of the bulk dataframe [{ds}] are not in the correct format! Skppping.')
            else:
                for i in ['Tapestri_format', 'condensed_format']:
                    if i in bulk_dfs[ds].columns:
                        bulk_dfs[ds].index = bulk_dfs[ds][i] # <--- set Tapestri_format as the index column
                    break
        print(f'[INFO] --- successfully loaded reference data {ds}')
        try:
            bulk_dfs[ds]['AF'] = bulk_dfs[ds]['t_alt_count'] / bulk_dfs[ds]['t_depth']
            bulk_dfs[ds] = bulk_dfs[ds][bulk_dfs[ds]['AF'] > 0] # only take vars with positive AF
            
        except KeyError:
            print(f'[INFO] --- AF cannot be calculated for {ds}')
            continue

    # -- b. matched bulk cohort 
    if bulk_info['matched_bulk_cohort'] is not None:
        if 'maf' not in bulk_info['matched_bulk_cohort'].keys():
            raise ValueError(f'[ERROR] --- for matched_bulk_cohort, {i} needs to be specified!')
        else:
            if 'tumor_sample_name' not in bulk_info['matched_bulk_cohort'].keys():
                print(f'[WARNING] --- for matched_bulk_cohort, tumor_sample_name is not specified. Not sample specific AF will be added.')
            bulk_dfs['matched_bulk_cohort'] = pd.read_csv(bulk_info['matched_bulk_cohort']['maf'], sep="\t")

            # get condensed variant format
            if not ('Tapestri_format' or 'condensed_format') in bulk_dfs['matched_bulk_cohort'].columns:
                try:
                    bulk_dfs['matched_bulk_cohort']['Tapestri_format'] = bulk_dfs['matched_bulk_cohort']['CHROM'] + ":" + bulk_dfs['matched_bulk_cohort']['POS'] + ":" + bulk_dfs['matched_bulk_cohort']['REF'] + "/" + bulk_dfs['matched_bulk_cohort']['ALT']
                    
                except KeyError:
                    print(f'[ERROR] --- the columns of the bulk dataframe [matched_bulk_cohort] are not in the correct format! Skppping.')
            bulk_dfs['matched_bulk_cohort'].index = bulk_dfs['matched_bulk_cohort']['Tapestri_format'] # <--- set Tapestri_format as the index column
            bulk_dfs['matched_bulk_cohort']['AF'] = bulk_dfs['matched_bulk_cohort']['t_alt_count'] / bulk_dfs['matched_bulk_cohort']['t_depth']
            bulk_dfs['matched_bulk_cohort'] = bulk_dfs['matched_bulk_cohort'][bulk_dfs['matched_bulk_cohort']['AF'] > 0] # only take vars with positive AF
            print(f'[INFO] --- successfully loaded [matched_bulk_cohort]. The DataFrame has {len(bulk_dfs["matched_bulk_cohort"].index.unique().tolist())} unique variants.')
            if bulk_info['matched_bulk_cohort']['annotation_summary_method'] is None:
                annotation_summary_method = 'max' # by default, the max vaf of all samples in cohort is annotated
            else:
                annotation_summary_method = bulk_info['matched_bulk_cohort']['annotation_summary_method']
    else:
        bulk_dfs['matched_bulk_cohort'] = None
        print(f'[INFO] --- matched_bulk_cohort is not given. Skipping.')

    # -- c. bulk normal        
    if bulk_info['matched_bulk_normal'] is not None:
        normal_sample_name = bulk_info['matched_bulk_normal']["normal_sample_name"]
        bulk_dfs['matched_bulk_normal'] = read_vcf_to_df(bulk_info['matched_bulk_normal']["vcf"], normal_sample_name)

        print(f'[INFO] --- successfully loaded [matched_bulk_normal]. The DataFrame has {len(bulk_dfs["matched_bulk_normal"].index.unique().tolist())} unique variants.')
        # AF calculation and index setting already in the `read_vcf_to_df` function
        # bulk_dfs['matched_bulk_normal']['AF'] = bulk_dfs['matched_bulk_normal']['t_alt_count'] / bulk_dfs['matched_bulk_normal']['t_depth']
    else:
        bulk_dfs['matched_bulk_normal'] = None
        print("[INFO] --- matched_bulk_normal is not given. Skipping.")

    # -- d. Tapestri PoN
    if bulk_info['tap_pon'] is not None:
        tap_pon = bulk_info['tap_pon']
        pon_dfs = {}
        for pon_v in tap_pon:
            pon_dfs[pon_v] = pd.read_csv(tap_pon[pon_v], sep="\t", index_col=0)
            # embed()
            if not check_matrix_format(pon_dfs[pon_v], CONDENSED_FORMAT):
                print(f'[ERROR] --- the index of {pon_v} is not in the correct condensed SNV format [chr]:[pos]:[ref]/[alt]. Skipping.')
                continue
            pon_dfs[pon_v]['median_sc_mut_prev'] = pon_dfs[pon_v].median(axis=1)
            print(f'[INFO] --- successfully loaded Tapestri PoN-{pon_v}. The DataFrame has {len(pon_dfs[pon_v].index.unique().tolist())} unique variants.')
    else:
        pon_dfs = None
        print('[INFO] --- tap_pon not given')

    # 2. ----- annotate VAR_DF -----

    # save the annotation stats to a dict
    annotation_stats = {}

    # helper functions
    def fetch_var_info(var: str, ref_df: pd.DataFrame, col_name: str) -> float:
        """
        fetch the AF of a variant from a reference dataframe
        
        Args:
            var: string, the variant in the format CHROM:POS:REF/ALT
            ref_df: pd.DataFrame, the reference dataframe
            col_name: string, the column info to extract
        Returns:
            AF: float, the AF of the variant
        """
        if var in ref_df.index:
            return ref_df.at[var, col_name]
        else:
            return np.nan
    fetch_var_info_v = np.vectorize(fetch_var_info, excluded = ['ref_df'])

    def compare_var_dfs(df_target, df_ref, df_ref_name, col_of_interest = None, summary_method = None) -> Tuple[pd.DataFrame, int, int]:
        """
        compare two dataframes and return the number of variants in each dataframe

        Args:
            df_target: pd.DataFrame, the target dataframe
            df_ref: pd.DataFrame, the reference dataframe
            df_ref_name: string, the name of the reference dataframe. This will be added to the output dataframe
            col_of_interest: string, the column of interest. If None, the presence/absence will be compared. Default: None
            summary_method: func, the summary method to use if multiple values are retrieved by df_ref.loc[var, col_of_interest]. Default: None

        Returns:
            # df_target_anntoated: pd.DataFrame, the target dataframe with the annotation
            ref_unique_var: int, the number of variants in the reference dataframe
            ref_unique_vars_covered: int, the number of variants in the reference dataframe covered by the target dataframe
        """
        ref_unique_count = len(df_ref.index.unique().tolist())
        ref_unique_count_covered = 0
        
        annotated_df = df_target.copy()
        for var_i in df_ref.index.unique().tolist():

            # get column of interest if present; else get presence/absence
            if col_of_interest is not None and col_of_interest in df_ref.columns:
                ann = df_ref.at[var_i, col_of_interest]
            else:
                print('[WARNING] --- column of interest not given/not present in df_ref column; using presence/absence.')
                ann = var_i in df_ref.index
            
            # if multiple records are present, summarize
            if type(ann) == pd.Series:
                if len(ann) > 1 and summary_method is None:
                    raise ValueError(f'[ERROR] --- variant {var_i} not unique in {ds}')
                if len(ann) > 1 and summary_method is not None:
                    ann = summary_method(ann)

            if not var_i in df_target.index:
                continue
            else:
                annotated_df.at[var_i, f'{col_of_interest}-{def_ref_name}'] = ann
                ref_unique_count_covered += 1

        return annotated_df, ref_unique_count, ref_unique_count_covered

    # (1). bulk data
    for ds in bulk_dfs:
        if bulk_dfs[ds] is None:
            continue

        # given the knowledge that the majority of vars in VAR_DF are artifacts, we can start from the bulk dfs instead of the VAR_DF to be annotated
        if ds == 'matched_bulk_cohort':
            # first, retrieve only the tumor sample
            if bulk_info['matched_bulk_cohort']['tumor_sample_name'] is not None:
                try:
                    sample_spec_df = bulk_dfs[ds].loc[bulk_dfs[ds]['Tumor_Sample_Barcode'] == bulk_info['matched_bulk_cohort']['tumor_sample_name']]
                    # embed()
                    VAR_DF, bulk_tumor_unique_vars_count, bulk_tumor_unique_vars_count_covered = compare_var_dfs(
                        VAR_DF, 
                        sample_spec_df, 
                        'matched_bulk_tumor'
                        'AF',
                        )
                    annotation_stats['bulk_tumor_unique_vars_count'] = bulk_tumor_unique_vars_count
                    annotation_stats['bulk_tumor_unique_vars_count_covered_by_Tapestri'] = bulk_tumor_unique_vars_count_covered
                    print(f"[INFO] --- finished annotating matched_bulk_tumor from matched_bulk_cohort")
                except KeyError:
                    print(f'[ERROR] --- Failed to retrieve {bulk_info["matched_bulk_cohort"]["tumor_sample_name"]}! Tumor sample-specific bulk AF not added')
            
            # second, get cohort-wide AF
            if annotation_summary_method == 'max':
                summary_method = np.max
            elif annotation_summary_method == 'min':
                summary_method = np.min
            elif annotation_summary_method == 'mean':
                summary_method = np.mean
            else:
                raise ValueError(f'[ERROR] --- Unknown annotation summary method: {annotation_summary_method}')

            VAR_DF, bulk_cohort_unique_vars_count, bulk_cohort_unique_vars_count_covered = compare_var_dfs(
                VAR_DF,
                bulk_dfs[ds],
                'matched_bulk_cohort',
                'AF',
                summary_method = summary_method
                )
            annotation_stats['bulk_cohort_unique_vars_count'] = bulk_cohort_unique_vars_count
            annotation_stats['bulk_cohort_unique_vars_count_covered_by_Tapestri'] = bulk_cohort_unique_vars_count_covered
            print(f"[INFO] --- finished annotating matched_bulk_cohort") 
        else:
            # each variant record should be unique in matched_bulk_tumor and broad_wes_pon
            VAR_DF, bulk_unique_var_count, bulk_unique_var_count_covered = compare_var_dfs(
                VAR_DF,
                bulk_dfs[ds],
                ds,
                'AF',
                )

            if ds == 'matched_bulk_tumor' or 'matched_bulk_normal':
                annotation_stats[f'{ds}_unique_vars_count'] = bulk_unique_var_count
                annotation_stats[f'{ds}_unique_vars_count_covered_by_Tapestri'] = bulk_unique_var_count_covered
            print(f"[INFO] --- finished annotating {ds}") 
            
    # (2). Tapestri PoN
    if pon_dfs is not None:
        for pon_vc in pon_dfs:
            # embed()
            # VAR_DF[f'median_sc_mut_prev-{pon_vc}'] = fetch_var_info_v(VAR_DF.index.values, ref_df = pon_dfs[pon_vc], col_name='median_sc_mut_prev')
            VAR_DF[f'median_sc_mut_prev-{pon_vc}'] = VAR_DF.index.map(lambda x: pon_dfs[pon_vc].at[x, 'median_sc_mut_prev'] if x in pon_dfs[pon_vc].index else np.nan)
            print(f"[INFO] --- finished annotating {pon_vc}")
            annotation_stats[f'{pon_vc}_unique_vars_count'] = len(pon_dfs[pon_vc].index.unique().tolist())
            annotation_stats[f'{pon_vc}_unique_vars_count_covered_by_Tapestri'] = int((VAR_DF[f'median_sc_mut_prev-{pon_vc}'] > 0).sum())
            # embed()

    # 3. ----- add info to H5 -----
    for col_name in VAR_DF.columns:
        sample_obj.dna.add_col_attr(col_name, VAR_DF[col_name].values)

    # 4. ----- write output -----
    mio.save(sample_obj, snakemake.output.annotated_h5)
    VAR_DF.to_csv(snakemake.output.tlod_metrics_tsv, sep="\t", index=True)
    with open(OUTPUT_DIR / f"{SAMPLE_NAME}-annotation_stats.json", "w") as out_f:
        json.dump(annotation_stats, out_f)