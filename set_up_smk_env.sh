# mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.32.3 --yes
# mamba activate snakemake
# # required packages on top of Snakemake: Plotly
# mamba install -c plotly plotly=5.13.1 --yes
# pip install -U kaleido # for plotly image export
# pip install pysam

# ----- set up LSF profile for Snakemake -----
# https://github.com/Snakemake-Profiles/lsf
# mamba install cookiecutter
# # create configuration directory that snakemake searches for profiles
# profile_dir="${HOME}/.config/snakemake"
# mkdir -p "$profile_dir"
# # use cookiecutter to create the profile in the config directory
# template="gh:Snakemake-Profiles/lsf"
# cookiecutter --output-dir "$profile_dir" "$template"

# set up the packages in resources/
WD=/data/iacobuzc/haochen/Tapestri_project/TapVarCallSmk # <------ edit this!!!
cd ${WD}/resources
# ----- custom packages -----
wget https://github.com/haochenz96/mosaic/archive/refs/tags/v2.0.0.zip -O mosaic-2.0.0.zip && unzip mosaic-2.0.0.zip && rm mosaic-2.0.0.zip
wget https://github.com/haochenz96/tea/archive/refs/tags/v2.0.zip -O tea-2.0.zip && unzip tea-2.0.zip && rm tea-2.0.zip
# append the above two packages to the conda env yaml file
echo "    - -e /data/iacobuzc/haochen/Tapestri_project/TapVarCallSmk/resources/mosaic-2.0.0" >> ${WD}/conda/envs/mosaic-custom.yaml
echo "    - -e /data/iacobuzc/haochen/Tapestri_project/TapVarCallSmk/resources/tea-2.0" >> ${WD}/conda/envs/mosaic-custom.yaml
# ----- BCFTOOLS -----
wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2 && tar -xvf bcftools-1.15.1.tar.bz2 && rm bcftools-1.15.1.tar.bz2
cd bcftools-1.15.1
make

# create necessary environments
cd ${WD}
snakemake \
	-s workflow/Snakefile \
	--profile lsf \
	--conda-prefix ${WD}/conda \
	--configfile ${WD}/configs/panel2_hg19-ucsc.yaml \
	--conda-create-envs-only 
