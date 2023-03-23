mamba create -c conda-forge -c bioconda -n snakemake snakemake

# required packages on top of Snakemake: Plotly
mamba install -c plotly plotly=5.13.1
pip install -U kaleido # for plotly image export