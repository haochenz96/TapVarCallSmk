conda activate gatk
gatk GetPileupSummaries \
    -I ${SC_BAM} \
    -V ${COMBINED_FILTERED_M2_BAM} \
    