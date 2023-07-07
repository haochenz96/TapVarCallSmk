rule step4_sc_mutect2_f_pass2:
    # scattered by single cell
    # this step, even at single-cell level, requires lots of resource. Might need some sample-specific hyperparameter tuning.
    input:
        sc_bam = "{sample_name}/1-sc_bams/{sample_name}_{cell_barcode}.bam",
        Q_VCF = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output: 
    	#vcf="$outpath/mutect2/__m2_sc_vcfs/${sample_name}_${cell_name}_somatic_m2.vcf.gz"
        sc_vcf = "{sample_name}/4-sc_mutect_f_call/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
        sc_vcf_idx = "{sample_name}/4-sc_mutect_f_call/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz.tbi",
        stats = "{sample_name}/4-sc_mutect_f_call/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz.stats", 
        # {output.stats} ensures Mutect2 finishes and is critical because this step usually takes long and several restarts
    params:
        REF = config['reference_info']['reference_genome'],
        mrpas = config['mutect2']['mrpas'],
        FILTERM2_OPS = parse_for_shell_args(config['mutect2']['filterm2_ops']),
    threads: 2
    resources: 
        mem_mb = lambda wildcards, attempt: ( 2**(attempt-1) ) * 2000,
        time_min = lambda wildcards, attempt: attempt * 179 + (attempt**2) * 60,
        # time_min = lambda wildcards, attempt: attempt * 159 + (attempt**2) * 60, # for second try
    conda:
        "../envs/mutect2.yaml"
    shell:
        """
    	sc_bai='{input.sc_bam}.bai'
    	if ! [ -f '$sc_bai' ]; then samtools index {input.sc_bam}; fi; 
        gatk Mutect2 \
            --java-options '-XX:-CreateCoredumpOnCrash' \
            -alleles {input.Q_VCF} \
            -L {input.Q_VCF} \
            --genotype-filtered-alleles \
            --max-reads-per-alignment-start {params.mrpas} \
            -R {params.REF} \
            -I {input.sc_bam} \
            -O {output.sc_vcf} && \
        samtools index {output.sc_vcf}
        """




rule step4_bcftools_pass2_filter_merge_sc_m2_f_calls:
    # scatter by sample
    # (1) for each single cell VCF, split multiallelic sites, apply only DP filter and index; 
    # (2) merge filtered single cell VCFs; no filter based on mutational prevalence across all cells; create statistics.
    input:
        sc_f_vcfs = get_step4_m2_f_vcfs,
    output:
        temp(directory("{sample_name}/4-sc_mutect_f_call/filtered_vcfs")),
        temp(directory("{sample_name}/4-sc_mutect_f_call/sc_TLOD_ann")),
        merged_f_vcf = "{sample_name}/4-sc_mutect_f_call/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        merged_f_vcf_stats = "{sample_name}/4-sc_mutect_f_call/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz.stats",
    params:
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'], # might want to use the latest release since the older version seems to have some bugs
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = 4000,
        time_min = lambda wildcards, attempt: attempt * 119,
    log: 
        "{sample_name}/logs/{sample_name}.snakemake_STEP4.log"
    shell:
        """
        echo "----------" >> {log}
        echo "[`date`]" >> {log}
        echo "[INFO] [STEP4] finished sc_mutect2_f_call for sample {wildcards.sample_name}" >> {log}

        # ----- (1) norm, filter, zip, index single-cell VCFs -----
        mkdir -p {wildcards.sample_name}/4-sc_mutect_f_call/filtered_vcfs
        for m2_vcf in {input.sc_f_vcfs}; do 
            
            out_name=$(basename ${{m2_vcf}})
            normed_vcf={wildcards.sample_name}/4-sc_mutect_f_call/filtered_vcfs/${{out_name/somatic_m2.vcf.gz/normed.vcf.gz}}

            {params.BCFTOOLS_EXEC} norm -m- ${{m2_vcf}} | \
            bgzip -c > \
            $normed_vcf

            tabix $normed_vcf
            echo cell -- ${{out_name}} -- norm, filter, zip, index done >> {log}

            ###############################################################################################    
            # @HZ 06/20/2022: 
            # ----- (2) we might want to keep the TLOD information emit by Mutect2 -----

            ANN_DIR={wildcards.sample_name}/4-sc_mutect_f_call/sc_TLOD_ann
            mkdir -p $ANN_DIR
            TLOD_ANN=$ANN_DIR/${{out_name/somatic_m2.vcf.gz/tlod.txt.gz}}
            
            # Extract INFO/TLOD into a tab-delimited annotation file
            {params.BCFTOOLS_EXEC} query \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%TLOD\n' \
                $normed_vcf | \
            bgzip -c \
                > $TLOD_ANN

            # Index the file with tabix
            tabix -s1 -b2 -e2 $TLOD_ANN

            # Create a header line for the new annotation
            if [ ! -f $ANN_DIR/header.txt ]; then
                echo -e '##FORMAT=<ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">' >> $ANN_DIR/header.txt
            fi

            # Transfer the annotation to sample and write new VCF
            smpl=$(zgrep -m1 ^#CHROM $normed_vcf | cut -f10)
            {params.BCFTOOLS_EXEC} annotate \
                -s $smpl \
                -a $TLOD_ANN \
                -h $ANN_DIR/header.txt \
                -c CHROM,POS,REF,ALT,FORMAT/TLOD \
                -Oz \
                -o ${{normed_vcf/normed.vcf.gz/normed.tlod.vcf.gz}} \
                $normed_vcf  
            tabix ${{normed_vcf/normed.vcf.gz/normed.tlod.vcf.gz}}
            ###############################################################################################
        done      

        # finally, merge
        {params.BCFTOOLS_EXEC} merge \
            --merge none \
            $(find {wildcards.sample_name}/4-sc_mutect_f_call/filtered_vcfs -name '*tlod.vcf.gz') | \
        {params.BCFTOOLS_EXEC} sort | \
        bgzip -c \
            > {output.merged_f_vcf}

        tabix {output.merged_f_vcf}

        {params.BCFTOOLS_EXEC} stats {output.merged_f_vcf} > {output.merged_f_vcf_stats} && \
         echo "[INFO] [STEP4] finished filtering and merging single-cell mutect2_f VCFs for sample {wildcards.sample_name}" >> {log}
        """

rule step4_bcftools_pass2_intersect_f_q_vcfs:
    # intersect q_vcf and f_vcf
    # all samples together
    input:
        f_vcf = "{sample_name}/4-sc_mutect_f_call/combined_vcf/{sample_name}-combined_VCF.filtered.vcf.gz",
        q_vcf = "{sample_name}/3-bcf_filter/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output:
        output_vcf = "{sample_name}/4-sc_mutect_f_call/{sample_name}-f_q_intersected.vcf.gz",
        output_vcf_stats = "{sample_name}/4-sc_mutect_f_call/{sample_name}-f_q_intersected.vcf.gz.stats",
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards, attempt: attempt * 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    log: 
        "{sample_name}/logs/{sample_name}.snakemake_STEP4.log"
    shell:
        """
        bcftools isec \
            -n=2 \
            -w2 \
            -Oz \
            -o {output.output_vcf} \
            {input.q_vcf} \
            {input.f_vcf}
        
        tabix {output.output_vcf}
        bcftools stats {output.output_vcf} > {output.output_vcf_stats} && \
         echo "[INFO] [STEP4] finished processing sample {wildcards.sample_name}" >> {log}
        """