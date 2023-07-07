'''
@HZ 07/04/2023: filter each single-sample H5, preparing for patient-wide genotyping:
- mininum prevalence of single cells called as mutated for an SNV
- maximum prevalence of single cells flagged as base_quality artifact for an SNV
- the number of cells called for an SNV should not be less than the number of cells flagged as base_quality artifact
'''

import argparse, sys, os
import logging
import mosaic.io as mio
from tea.format import CONDENSED_SNV_FORMAT, check_matrix_format
import pandas as pd

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

def main(args):

    # ----- read input H5 file and bq flag file -----
    sample_obj = mio.load(args.input_h5)
    sample_obj.dna.genotype_variants(
        het_vaf=20, 
        hom_vaf=80, 
        min_dp=8, 
        min_alt_read = 3, 
        min_gq = -1, 
        assign_low_conf_genotype=False
    )
    num_cells = sample_obj.dna.shape[0]
    num_snvs = sample_obj.dna.shape[1]
    logging.info(f"[STEP5-prev_filter] Number of cells: {num_cells}")
    logging.info(f"[STEP5-prev_filter] Number of SNVs: {num_snvs}")

    mut_prev_df = sample_obj.dna.get_attribute('mut_filtered', constraint='row').sum(axis=0)
    bq_prev_df = pd.read_csv(args.bq_info_csv, index_col=0, header=None, names=['bq_sc_prev'])
    if not check_matrix_format(bq_prev_df, CONDENSED_SNV_FORMAT):
        exit(1, "ERROR: bq_info_csv is not in the correct format.")
    
    mut_prev_df.name='mut_sc_prev'
    snv_prev_df = mut_prev_df.to_frame()
    snv_prev_df['bq_sc_prev'] = snv_prev_df.index.map(bq_prev_df['bq_sc_prev']).fillna(0)

    # ----- filter -----
    # get params:
    mut_prev_threshold = args.mut_prev_threshold
    bq_prev_threshold = args.bq_prev_threshold

    # 1. mutational prevalence
    mask = snv_prev_df['mut_sc_prev'] >= (num_cells * mut_prev_threshold)
    print(mask.sum())
    logging.info(f"[STEP5-prev_filter] filtered SNVs with mutational prevalence < {mut_prev_threshold}: {(snv_prev_df['mut_sc_prev'] < (num_cells*mut_prev_threshold)).sum()}")

    # 2. bq prevalence
    mask = mask & (snv_prev_df['bq_sc_prev'] < (num_cells * bq_prev_threshold))
    print(mask.sum())
    logging.info(f"[STEP5-prev_filter] filtered SNVs flagged as base_quality artifact with prevalence >= {bq_prev_threshold}: {(snv_prev_df['bq_sc_prev'] >= num_cells * bq_prev_threshold).sum()}")
    
    # 3. num of cells called >= num of cells flagged as bq artifact
    mask = mask & (snv_prev_df['mut_sc_prev'] > snv_prev_df['bq_sc_prev'])
    print(mask.sum())
    logging.info(f"[STEP5-prev_filter] filtered SNVs flagged as base_quality artifact in more cells than they are confidently called: {(snv_prev_df['mut_sc_prev'] <= snv_prev_df['bq_sc_prev']).sum()}")

    voi = snv_prev_df[mask].index.tolist()
    logging.info(f"[STEP5-prev_filter] Number of SNVs left after filtering: {len(voi)}")

    # ----- write output list file -----
    if args.output_list_f is not None:
        pd.DataFrame(index=voi).to_csv(args.output_list_f, header=False, index=True)
        logging.info(f"[STEP5-prev_filter] Written to output list file: {args.output_list_f}")

    # ----- write output H5 -----
    if args.output_h5 is not None:
        sample_obj.dna = sample_obj.dna[sample_obj.dna.barcodes(), voi]
        mio.save(sample_obj, args.output_h5)
        logging.info(f"[STEP5-prev_filter] Written to output H5 file: {args.output_h5}")
    # else:
    #     logging.info(f"[STEP5-prev_filter] No output H5 file is written.")

    # ----- construct output VCF file -----
    if args.output_vcf is not None:
        df = pd.DataFrame(index = voi)

        df['#CHROM'] = df.index.str.split(':').str[0]
        # df['chr_num'] = df['#CHROM'].str.strip('chr').astype(int)

        df['POS'] = df.index.str.split(':').str[1].astype(int)
        # df.sort_values(['chr_num','POS'], inplace=True)
        df['ID'] = '.'
        df['REF'] = df.index.str.split(':').str[2].str.split('/').str[0]
        df['ALT'] = df.index.str.split(':').str[2].str.split('/').str[1]
        df['QUAL'] = '.'
        df['FILTER'] = '.'
        df['INFO'] = '.'

        with open(args.output_vcf, 'w') as out:
            out.write('##fileformat=VCFv4.2\n')
            out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

            for _, row in df.iterrows():
                out.write(f"{row['#CHROM']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\n")
            
            # elif args.genome == 'hg19-b37':
            #     for _, row in df.iterrows():
            #         out.write(f"{row['chr_num']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\n")
        logging.info(f"[STEP5-prev_filter] Output VCF file written to {args.output_vcf}")
    # --------------------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_h5', type=str, help='path to write output H5.', required=True)
    parser.add_argument('--bq_info_csv', type=str, help='base-quality info output by the custom pipeline.', required=True)
    parser.add_argument('--output_list_f', type=str, help='path to write output list of SNVs.', default=None)
    parser.add_argument('--output_vcf', type=str, help='path to write output vcf.', default=None)
    parser.add_argument('--output_h5', type=str, help='path to write output vcf.', default=None)
    parser.add_argument('--mut_prev_threshold', type=float, help='minimum mutation prevalence threshold.', default=0.005)
    parser.add_argument('--bq_prev_threshold', type=float, help='maximum prevalence of cells flagged as base quality artifact.', default=0.005)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)