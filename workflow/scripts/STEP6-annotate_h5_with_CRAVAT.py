# @HZ 06/20/2022
# script to add SNV annotation to H5 file

# first, install customized-mosaic and tea from Github

from pathlib import Path
import mosaic.io as mio
from tea.utils import *
from tea.plots import *
from tea.cravat import *
import requests
import base64
import json

# ----- fetch inputs -----
H5 = snakemake.input.bulk_annotated_h5
sample_name = snakemake.wildcards.sample_name
cravat_annotators = snakemake.params.cravat_annotators

with open (snakemake.params.log_file, 'a') as f:
    sys.stdout = f

    # 1. establish connection to OpenCRAVAT server:
    new_session = requests.Session()

    # log in with credentials 
    reply = new_session.get(
        'https://run.opencravat.org/server/login', 
        headers={'Authorization': 'Basic ' + base64.b64encode(b'zhangh5@mskcc.org:Hfhd68MdTtHn5UE').decode()})

    if not reply.json() == 'success':
        print('[Warning] CRAVAT login failed!')
        exit(1)
    
    sample_obj = mio.load(H5, sample_name = sample_name, raw = False)
    sample_obj.dna.genotype_variants(
        min_dp = 5,
        min_alt_read = 3
    )
    # mut_prev_array = sample_obj.dna.get_attribute(
    #     'mut', constraint='row'
    # ).sum(axis=0)
    # vars_of_interest = mut_prev_array.index[mut_prev_array >= 10] # <----- take any variant present in more than 10 cells

    # CRAVAT input will be written to "wd / CRAVAT"
    (wd / 'CRAVAT' / 'CRAVAT_input').mkdir(exist_ok=True, parents=True)
    (wd / 'CRAVAT' / 'CRAVAT_output').mkdir(exist_ok=True, parents=True)
    cravat_input = wd / 'CRAVAT' / 'CRAVAT_input' / f'{sample_i}_CRAVAT_input.txt'
    cravat_output = wd / 'CRAVAT' / 'CRAVAT_output' / f'{sample_i}_CRAVAT_output.txt'
    write_cravat_input(vars_of_interest, cravat_input, sample_i, mut_prev_array)
    if not cravat_output.is_file():
        cravat_params = json.dumps({
            "annotators": cravat_annotators, 
            "reports": ["text"], 
            "assembly": "hg19", 
            "note": f"{sample_name}-analysis",
            })

        post = new_session.post(
            'https://run.opencravat.org/submit/submit', 
            files={'file_0': open(cravat_input)}, 
            data={'options': cravat_params
                }
        )
        
        cravat_job_id = post.json()['id']
        get_cravat_output(new_session, cravat_job_id, cravat_output)
        print(f'[INFO] Finished CRAVAT analysis for sample -- {sample_i}.')
    else:
        print('[WARNING] CRAVAT output file already exists! Skipping...')
