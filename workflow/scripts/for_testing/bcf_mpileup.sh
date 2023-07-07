conda activate samtools

# ----- get variables -----
# inputs
SC_BAM=/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/RA16_29/RA16_29-16_1/1-sc_bams/RA16_29-16_1_cell_221.bam
CANDIDATE_ALLELE=/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/RA16_29/RA16_29-16_1/4-bcf_genotyping/RA16_29-16_1-candidate_alleles.tsv.gz

# outputs
OUT_DIR=/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/RA16_29/test
SC_MPILEUP_VCF=${OUT_DIR}/test.vcf.gz

# params
REF_GENOME=/home/zhangh5/work/commons/Tapestri_ref/hg19-b37/hg19-b37.fasta
PANEL_AMPLICON=/home/zhangh5/work/Tapestri_batch2/panel_3359/panel_3359_hg19-b37/3359.amplicons.bed

# # create a candidate allele file
# bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $Q_VCF | \
# bgzip -c > $CANDIDATE_ALLELE && \
# tabix -s1 -b2 -e2 $CANDIDATE_ALLELE

# mpileup and call
bcftools mpileup \
    -Ou \
    -R ${PANEL_AMPLICON} \
    -f ${REF_GENOME} \
    --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
    --max-depth 100000 \
    --max-idepth 100000 \
    --no-BAQ \
    ${SC_BAM} | \
bcftools call \
    --keep-alts \
    -C alleles \
    -T ${CANDIDATE_ALLELE} \
    --multiallelic-caller \
    -Oz \
    -o ${SC_MPILEUP_VCF} && \
tabix ${SC_MPILEUP_VCF}


# extract INFO/AD into a tab-delimited annotation file
SAMPLE_NAME=$(bcftools query -l $SC_MPILEUP_VCF)

# only get the first two alleles
bcftools query -f '%CHROM\t%POS\t%AD{0},%AD{1}\n' ${SC_MPILEUP_VCF} | bgzip -c > ${OUT_DIR}/annot.txt.gz && \
tabix -s1 -b2 -e2 ${OUT_DIR}/annot.txt.gz

echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Total allelic depths (high-quality bases)">' >> ${OUT_DIR}/hdr.txt
bcftools annotate \
    -s ${SAMPLE_NAME} \
    -a ${OUT_DIR}/annot.txt.gz \
    -h ${OUT_DIR}/hdr.txt \
    -c CHROM,POS,FORMAT/AD \
    -Oz \
    -o ${OUT_DIR}/test.ann.vcf.gz \
    ${SC_MPILEUP_VCF}

# # FOR TESTING ONLY
# bcftools mpileup \
#     -Ou \
#     -R ${PANEL_AMPLICON} \
#     -f ${REF_GENOME} \
#     --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
#     --max-depth 100000 \
#     --max-idepth 100000 \
#     --no-BAQ \
#     -Oz \
#     -o test.mpileup.vcf.gz \
#     ${SC_BAM} 