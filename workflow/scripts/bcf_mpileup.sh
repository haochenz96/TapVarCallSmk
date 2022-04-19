# ----- get variables -----
# inputs
SC_BAM=snakemake.input.SC_BAM
CANDIDATE_ALLELE=snakemake.input.CANDIDATE_ALLELE

# outputs
SC_MPILEUP_VCF=snakemake.output.SC_MPILEUP_VCF

# params
REF_GENOME=snakemake.params.REF_GENOME
PANEL_AMPLICON=snakemake.params.PANEL_AMPLICON

# # create a candidate allele file
# bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $Q_VCF | \
# bgzip -c > $CANDIDATE_ALLELE && \
# tabix -s1 -b2 -e2 $CANDIDATE_ALLELE

# mpileup and call
bcftools mpileup \
    -Ou \
    -R $PANEL_AMPLICON \
    -f $REF_GENOME \
    --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
    -d 100000 \
    $SC_BAM | \
bcftools call \
    --keep-alts \
    -C alleles \
    -T $CANDIDATE_ALLELE \
    --multiallelic-caller \
    -Ov - > $SC_MPILEUP_VCF && \
bgzip $SC_MPILEUP_VCF && \
tabix ${SC_MPILEUP_VCF}.gz