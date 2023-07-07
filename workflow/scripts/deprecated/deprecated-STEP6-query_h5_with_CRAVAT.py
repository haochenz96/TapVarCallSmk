# @HZ 06/20/2022
# script to query variants from a PoN_v1_filtered H5 with CRAVAT
# inputs:
#   - an H5 file
# outputs:
#   - a file with the variants and their CRAVAT annotations

# first, install customized-mosaic and tea from Github

from pathlib import Path
import mosaic.io as mio
from tea.utils import *
from tea.plots import *
from tea.cravat import *
import requests
import base64
import json
import sys, argparse, yaml

# for testing
# SAMPLE_NAME = 'M13-1'
# H5 = '/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/M13-1_combined_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5'
# __output_dir = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/annotations')
# CRAVAT_OUTPUT = __output_dir / 'M13-1_ponv1_filtered_CRAVAT_output.txt'
# LOG_FILE = __output_dir / 'STEP6_CRAVAT_log.txt'
# OPENCRAVAT_USERNAME = 'zhangh5@mskcc.org'
# OPENCRAVAT_PASSWORD = 'Hfhd68MdTtHn5UE'
# CRAVAT_ANNOTATORS = ["clinvar", "gnomad3", "chasmplus", "chasmplus_PAAD","civic","cosmic","dbsnp","dbsnp_common","clinvar","gnomad3",'thousandgenomes','go','ndex']
# CRAVAT_INPUT_PATH = __output_dir / 'M13-1_ponv1_filtered_CRAVAT_input.txt'

def main(args):
    # ----- io -----
    SAMPLE_NAME = args.sample_name
    # --- inputs
    H5 = args.input_h5
    # --- outputs
    if args.output_dir is None:
        CRAVAT_OUTPUT = args.cravat_output
        CRAVAT_INPUT_PATH = CRAVAT_OUTPUT.parent / f'{SAMPLE_NAME}_CRAVAT_input.txt'
    else:
        __output_dir = Path(args.output_dir)
        __output_dir.mkdir(exist_ok=True, parents=True)
        CRAVAT_OUTPUT = __output_dir / f'{SAMPLE_NAME}_CRAVAT_output.txt'
        CRAVAT_INPUT_PATH = __output_dir / f'{SAMPLE_NAME}_CRAVAT_input.txt'
    # --- params
    ###############################
    LOG_FILE = args.log_file
    if LOG_FILE is None:
        logging.basicConfig(
            stream=sys.stdout,
            level=logging.INFO, 
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    else:
        logging.basicConfig(
            filename = LOG_FILE,
            level = logging.INFO,
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    ###############################  
    cravat_settings_yaml = args.cravat_settings_yaml
    with open(cravat_settings_yaml, "r") as f:
        cravat_settings = yaml.safe_load(f)

    try:
        OPENCRAVAT_USERNAME = cravat_settings['username']
        OPENCRAVAT_PASSWORD = cravat_settings['password']
        CRAVAT_ANNOTATORS = cravat_settings['annotators']
    except KeyError:
        logging.error('Please check your cravat_settings.yaml file')
        sys.exit(1)
    
    # 1. establish connection to OpenCRAVAT server:
    new_session = requests.Session()

    # log in with credentials 
    cred_str = base64.b64encode(bytes(OPENCRAVAT_USERNAME + ":" + OPENCRAVAT_PASSWORD, 'utf-8')).decode('utf-8')
    reply = new_session.get(
        'https://run.opencravat.org/server/login', 
        headers={'Authorization': 'Basic ' + cred_str})
    # original cred_str ----->> base64.b64encode(b'zhangh5@mskcc.org:Hfhd68MdTtHn5UE').decode() 

    if not reply.json() == 'success':
        logging.warning('CRAVAT login failed!')
        exit(1)
    
    sample_obj = mio.load(H5, name = SAMPLE_NAME, raw = False)
    sample_obj.dna.genotype_variants(
        min_dp = 8,
        min_alt_read = 3
    )
    mut_prev_array = sample_obj.dna.get_attribute(
        'mut_filtered', constraint='row'
    ).sum(axis=0)
    # vars_of_interest = sample_obj.dna.ids()[mut_prev_array > 0.01 * sample_obj.dna.shape[0]] # <----- take any variant present in more than 1% cells
    # vars_of_interest = sample_obj.dna.ids()
    vars_of_interest = sample_obj.dna.ids()[
            (np.isnan(sample_obj.dna.col_attrs['mean_sc_mut_prev-tap_pon_v2']) | 
            ~np.isnan(sample_obj.dna.col_attrs['AF-matched_bulk_normal'])) & (mut_prev_array >= 3)
        ] # filter out ponv2 variants as well
        
    write_cravat_input(vars_of_interest, CRAVAT_INPUT_PATH, SAMPLE_NAME, mut_prev_array)
    logging.info(f'{len(vars_of_interest)} variants of interest written to CRAVAT input')
    if not CRAVAT_OUTPUT.is_file():
        cravat_params = json.dumps({
            "annotators": CRAVAT_ANNOTATORS, 
            "reports": ["text"], 
            "assembly": "hg19", 
            "note": f"{SAMPLE_NAME}-analysis",
            })

        post = new_session.post(
            'https://run.opencravat.org/submit/submit', 
            files={'file_0': open(CRAVAT_INPUT_PATH)}, 
            data={'options': cravat_params
                }
        )
        
        cravat_job_id = post.json()['id']
        get_cravat_output(new_session, cravat_job_id, CRAVAT_OUTPUT)
        logging.info(f'Finished CRAVAT analysis for sample -- {SAMPLE_NAME}.')
    else:
        logging.warning('CRAVAT output file already exists! Skipping...')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type = str, help='sample name')
    parser.add_argument('--input_h5', type = str, help='input_h5')
    parser.add_argument('--cravat_settings_yaml', type = str, required=True, help='bulk_info_yaml')
    parser.add_argument('--cravat_output', type = str, help='output CRAVAT report')
    parser.add_argument('--output_dir', type = str, help='output_dir; overrides the output file name if specified.')
    parser.add_argument('--log_file', type = str, default=None, help='log file')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)