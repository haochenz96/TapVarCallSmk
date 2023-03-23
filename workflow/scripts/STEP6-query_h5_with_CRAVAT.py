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

# ----- io -----
SAMPLE_NAME = snakemake.wildcards.sample_name
# --- inputs
H5 = snakemake.input.ponv1_filtered_h5
# --- outputs
CRAVAT_OUTPUT = Path(snakemake.output.cravat_output)
# --- params
LOG_FILE = Path(snakemake.params.log_file)
OPENCRAVAT_USERNAME = snakemake.params.opencravat_username
OPENCRAVAT_PASSWORD = snakemake.params.opencravat_password
CRAVAT_ANNOTATORS = snakemake.params.cravat_annotators
CRAVAT_INPUT_PATH = Path(snakemake.params.cravat_input_path)

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

with open (LOG_FILE, 'a') as f:
    sys.stdout = f # redirect stdout to log file

    # 1. establish connection to OpenCRAVAT server:
    new_session = requests.Session()

    # log in with credentials 
    cred_str = base64.b64encode(bytes(OPENCRAVAT_USERNAME + ":" + OPENCRAVAT_PASSWORD, 'utf-8')).decode('utf-8')
    reply = new_session.get(
        'https://run.opencravat.org/server/login', 
        headers={'Authorization': 'Basic ' + cred_str})
    # original cred_str ----->> base64.b64encode(b'zhangh5@mskcc.org:Hfhd68MdTtHn5UE').decode() 

    if not reply.json() == 'success':
        print('[Warning] CRAVAT login failed!')
        exit(1)
    
    sample_obj = mio.load(H5, name = SAMPLE_NAME, raw = False)
    sample_obj.dna.genotype_variants(
        min_dp = 5,
        min_alt_read = 3
    )
    mut_prev_array = sample_obj.dna.get_attribute(
        'mut', constraint='row'
    ).sum(axis=0)
    vars_of_interest = mut_prev_array.index # <----- take any variant present in more than 10 cells
    
    write_cravat_input(vars_of_interest, CRAVAT_INPUT_PATH, SAMPLE_NAME, mut_prev_array)
    print(f'[Info] {len(vars_of_interest)} variants of interest written to CRAVAT input')
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
        print(f'[INFO] Finished CRAVAT analysis for sample -- {SAMPLE_NAME}.')
    else:
        print('[WARNING] CRAVAT output file already exists! Skipping...')
