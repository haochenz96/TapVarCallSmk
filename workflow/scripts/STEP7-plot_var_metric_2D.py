# @HZ 06/27/2022
# script to plot 2D scatterplot of of TLOD metrics of Tapestri data
# do this in the mosaic-custom conda environment
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
import argparse
import yaml
from pathlib import Path
from IPython import embed

# ----- io ----- 
H5 = snakemake.input.ponv1_filtered_h5
SAMPLE_NAME = snakemake.wildcards.sample_name
FIG1 = snakemake.output.fig1
FIG2 = snakemake.output.fig2
FIG3 = snakemake.output.fig3
OUTPUT_DIR = Path(snakemake.params.output_dir)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# ----- for testing -----
# H5 = '/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/RA16_29/RA16_29-32_1/OUTPUTS_from_m2_f/RA16_29-32_1_DNA_CNV_m2_f.cravat_annotated.h5'
# SAMPLE_NAME = 'RA16_29-32_1'
# OUTPUT_DIR = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/RA16_29/RA16_29-32_1/OUTPUTS_from_m2_f/test')

sample_obj = mio.load(H5, raw=False)
sample_obj.dna.genotype_variants(
    min_dp = 5,
    min_alt_read = 3,
)
mut_prev_array = sample_obj.dna.get_attribute(
    'mut', constraint='row'
).sum(axis=0)
TLOD_METRIC_DF = pd.DataFrame(index = mut_prev_array.index)
TLOD_METRIC_DF['sc_mut_prev'] = mut_prev_array.values
for col in sample_obj.dna.col_attrs.keys():
    if not any(kwd in col for kwd in ['bulk', 'pon','TLOD']):
        continue
    else:
        print(f'[INFO] adding {col}')
        TLOD_METRIC_DF[col] = sample_obj.dna.col_attrs[col]

# subsample = args.subsample
# if subsample is not None:
#     TLOD_METRIC_DF = TLOD_METRIC_DF.sample(frac=float(subsample), axis=0) # randomly subsample for efficiency

bulk_data_cols = [i for i in TLOD_METRIC_DF.columns if 'bulk' in i]
print(f'[INFO] bulk_data_cols: {bulk_data_cols}')
TLOD_METRIC_DF['bulk_validation'] = (TLOD_METRIC_DF[bulk_data_cols]>0).any(axis=1)
TLOD_METRIC_DF['bulk_validation'] = TLOD_METRIC_DF['bulk_validation'].map(
    {True: 'Yes', False: 'No'}
)

TLOD_METRIC_DF['ponv2_presence'] = (TLOD_METRIC_DF['median_sc_mut_prev-tap_pon_v2']>0).map(
    {True: 'Yes', False: 'No'}
)

fig = px.scatter(
    TLOD_METRIC_DF,
    # x = 'sc_mean_TLOD',
    x = 'sc_mut_prev',
    y = 'sc_max_TLOD',
    color = 'bulk_validation',
    marginal_y = 'violin',
    color_discrete_map={
        'Yes': "green",
        'No': "grey",
        },
    log_x = True,
    log_y = True,
    opacity = 0.6,
    size_max = 2,
)

# fig.update_xaxes(
#     title = f'%amplicons with mean_rc >= {rc_thres}',
# )
# fig.update_layout(
#     title=f"""
#     {sample_name} single-cell amplicon performance distribution<br>
#     <sup>candidates: {total_cell_num} cells</sup>
#     <sup>cells called: {ncalled_cells} cells</sup>
#     """,
# )
# fig.add_vline(
#     80, 
#     line_color='red', 
#     line_width=2, row=1, col=1, opacity=0.2,
#     annotation_text=f"default cellfinder cutoff = 80%",
#     annotation_font_size=10,
#     annotation_font_color="red",
# )

fig.write_image(FIG1, width=500, height=500, scale=2)
print(f'[SUCCESS] finished plotting 2D scatterplot of TLOD metrics of Tapestri data, colored by bulk validation')

fig2 = px.scatter(
    TLOD_METRIC_DF,
    # x = 'sc_mean_TLOD',
    x = 'sc_mut_prev',
    y = 'sc_max_TLOD',
    color = 'ponv2_presence',
    marginal_y = 'violin',
    color_discrete_map={
        'No': "grey",
        'Yes': "red",
        },
    log_x = True,
    log_y = True,
    opacity = 0.6,
    size_max = 2,
)

fig2.write_image(FIG2, width=500, height=500, scale=2)
print(f'[SUCCESS] finished plotting 2D scatterplot of TLOD metrics of Tapestri data, colored by PoN V2 presence')

fig3 = px.scatter(
    TLOD_METRIC_DF,
    # x = 'sc_mean_TLOD',
    x = 'sc_mut_prev',
    y = 'sc_mean_TLOD',
    color = 'ponv2_presence',
    marginal_y = 'violin',
    color_discrete_map={
        'No': "grey",
        'Yes': "red",
        },
    log_x = True,
    log_y = True,
    opacity = 0.6,
    size_max = 2,
)
fig3.write_image(FIG3, width=500, height=500, scale=2)
print(f'[SUCCESS] finished plotting 2D scatterplot of TLOD metrics of Tapestri data, colored by PoN V2 presence')