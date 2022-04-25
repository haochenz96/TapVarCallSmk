rule generate_candidate_allele:
    # generate a candidate allele from Q_VCF
    input:
        Q_VCF = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output: 
    	CANDIDATE_ALLELE = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.tsv.gz",
    threads: lambda wildcards, attempt: attempt * 2
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {input.Q_VCF} | "
        "bgzip -c > {output.CANDIDATE_ALLELE} && "
        "tabix -s1 -b2 -e2 {output.CANDIDATE_ALLELE} "

rule sc_mpileup:
    # scattered by single cell
    # mpileup and call variants 
    input:
        SC_BAM = "{sample_name}/1-sc_bams/{sample_name}_{cell_num_index}.bam",
        CANDIDATE_ALLELE = "{sample_name}/4-bcf_genotyping/{sample_name}-candidate_alleles.tsv.gz",
    output:
        SC_MPILEUP_VCF = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}/{sample_name}_{cell_num_index}_mpileup.vcf.gz",
        SC_MPILEUP_VCF_raw = "{sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{sample_name}_{cell_num_index}/{sample_name}_{cell_num_index}_raw_counts.vcf.gz",
    params:
        REF_GENOME = config['reference_info']['reference_genome'],
        PANEL_AMPLICON = config['reference_info']['panel_amplicon_file'],
    threads: lambda wildcards, attempt: attempt * 2
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 119,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        OUT_DIR="{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{wildcards.sample_name}_{wildcards.cell_num_index}"
        mkdir -p $OUT_DIR

        # @STEP1 mpileup and call
        # @HZ 04/24/2022 added option "--no-BAQ"; try to remove base quality consideration to recover MNP's
        bcftools mpileup \
            -Ou \
            -R {params.PANEL_AMPLICON} \
            -f {params.REF_GENOME} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            --max-depth 100000 \
            --max-idepth 100000 \
            --no-BAQ \
            {input.SC_BAM} | \
        bcftools call \
            --keep-alts \
            -C alleles \
            -T {input.CANDIDATE_ALLELE} \
            --multiallelic-caller \
            -Oz \
            -o {output.SC_MPILEUP_VCF} && \
        tabix {output.SC_MPILEUP_VCF}

        # @STEP2 get raw counts
        # extract INFO/AD into a tab-delimited annotation file
        SAMPLE_NAME=$(bcftools query -l {output.SC_MPILEUP_VCF})

        bcftools query -f '%CHROM\t%POS\t%AD{{0}},%AD{{1}}\n' {output.SC_MPILEUP_VCF} | bgzip -c > $OUT_DIR/annot.txt.gz && \
        tabix -s1 -b2 -e2 $OUT_DIR/annot.txt.gz

        echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Total allelic depths (high-quality bases)">' >> $OUT_DIR/hdr.txt
        bcftools annotate \
            -s $SAMPLE_NAME \
            -a $OUT_DIR/annot.txt.gz \
            -h $OUT_DIR/hdr.txt \
            -c CHROM,POS,FORMAT/AD \
            -Oz \
            -o {output.SC_MPILEUP_VCF_raw} \
            {output.SC_MPILEUP_VCF}

        # @STEP3 calculate AF
        """

rule filter_and_merge_sc_mpileup:
    # scatter by sample
    # (1) for each single cell VCF, split multiallelic sites, apply only DP filter and index; 
    # (2) merge filtered single cell VCFs; no filter based on mutational prevalence across all cells; create statistics.
    input:
        sc_mpileup_vcfs = get_step4_mpileup_vcfs,
    output:
        merged_mpileup_vcf = "{sample_name}/4-bcf_genotyping/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        merged_mpileup_vcf_stats = "{sample_name}/4-bcf_genotyping/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz.stats",
    params:
        sc_vcf_filters = parse_bcf_filters(config['bcftools']['pass2_sc_vcf_filters']),
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    shell:
        """
        mkdir -p {wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_filtered_vcfs
        for m2_vcf in {input.sc_mpileup_vcfs}; do 
            
            out_name=$(basename ${{m2_vcf}})
            
            
            bcftools norm -m- ${{m2_vcf}} | \
            bcftools filter -i '{params.sc_vcf_filters}' | \
            bgzip -c > \
            {wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_filtered_vcfs/${{out_name/.vcf.gz/.hard_filtered.vcf.gz}}

            tabix {wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_filtered_vcfs/${{out_name/.vcf.gz/.hard_filtered.vcf.gz}}
            echo cell -- ${{out_name}} -- norm, filter, zip, index done
        done

        bcftools merge --merge none $(find {wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_filtered_vcfs -name '*hard_filtered.vcf.gz') | \
        bcftools sort | \
        bgzip -c > \
        {output.merged_mpileup_vcf}

        tabix {output.merged_mpileup_vcf}

        bcftools stats {output.merged_mpileup_vcf} > {output.merged_mpileup_vcf_stats}

        """