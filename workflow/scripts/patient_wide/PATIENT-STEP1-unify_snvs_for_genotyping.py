'''
@HZ 07/05/2023: given several lists of SNVs, union them and format as a VCF to be genotyped
- mininum prevalence of single cells called as mutated for an SNV
- maximum prevalence of single cells flagged as base_quality artifact for an SNV
- the number of cells called for an SNV should not be less than the number of cells flagged as base_quality artifact
'''

import argparse, sys, os
import logging
import pandas as pd

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

def main(args):

    # ----- read in SNV lists -----
    dfs = [pd.read_csv(f, index_col=0) for f in args.input_snv_lists]
    logging.info(f"[STEP1-unify_SNVs] {len(dfs)} SNV lists read in.")
    snv_sets = [set(df.index) for df in dfs]
    voi = list(snv_sets[0].union(*snv_sets[1:]))
    logging.info(f"[STEP1-unify_SNVs] {len(voi)} unique SNVs in total.")

    # ----- construct output VCF file -----
    output_df = pd.DataFrame(index = voi)
    output_df['#CHROM'] = output_df.index.str.split(':').str[0]
    output_df['chr_num'] = output_df['#CHROM'].str.strip('chr')

    # custom key method for sorting, per https://stackoverflow.com/questions/72801110/sorting-a-list-of-chromosomes-in-the-correct-order
    key={s:i for i,s in enumerate([str(x) for x in list(range(1,22))+['X','Y']],1)}
    output_df['chr_num'] = output_df['chr_num'].map(key)

    output_df['POS'] = output_df.index.str.split(':').str[1].astype(int)
    output_df.sort_values(['chr_num','POS'], inplace=True)
    output_df['ID'] = '.'
    output_df['REF'] = output_df.index.str.split(':').str[2].str.split('/').str[0]
    output_df['ALT'] = output_df.index.str.split(':').str[2].str.split('/').str[1]
    output_df['QUAL'] = '.'
    output_df['FILTER'] = '.'
    output_df['INFO'] = '.'

    with open(args.output_vcf, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for _, row in output_df.iterrows():
            out.write(f"{row['#CHROM']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\n")
        
        # elif args.genome == 'hg19-b37':
        #     for _, row in output_df.iterrows():
        #         out.write(f"{row['chr_num']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\n")
    logging.info(f"[STEP5-prev_filter] Output VCF file written to {args.output_vcf}")
    # --------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_snv_lists', nargs="*", help='path to input SNV s.', required=True)
    parser.add_argument('--output_vcf', type=str, help='path to write output vcf.', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)