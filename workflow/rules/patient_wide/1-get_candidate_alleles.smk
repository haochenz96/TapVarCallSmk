rule patient_step1_gather_individual_snv_lists:
    input:
        sample_snv_lists = expand(
            "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.mpileup.prev_filtered.csv", sample_name = config['patient_info']['sample_names']
            ),
    output:
        patient_union_vcf = "fillout/input/{patient_name}.snv_union.vcf",
    params:
        script = os.path.join(scripts_dir, "patient_wide/PATIENT-STEP1-unify_snvs_for_genotyping.py")
    log:
        "fillout/logs/{patient_name}-patient-wide_genotype.log"
    shell:
        """
        python {params.script} \
            --input_snv_lists {input.sample_snv_lists} \
            --output_vcf {output.patient_union_vcf} \
            > {log} 2>&1
        """

rule patient_step1_generate_candidate_alleles_for_bcf:
    # generate a candidate allele from Q_VCF
    # @HZ 08/28/2022: need to merge multiallelic sites for the next step (mpileup+call) to work
    input:
        patient_union_vcf = "fillout/input/{patient_name}.snv_union.vcf",
    output: 
        CANDIDATE_ALLELE = "fillout/input/{patient_name}.snv_union.for_genotyping.vcf.gz",
        CANDIDATE_ALLELE_for_py = "fillout/input/{patient_name}.snv_union.multiallelic.for_py.csv",
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../../envs/bcftools.yaml",
    params:
        Q_VCF_atomized = "fillout/input/{patient_name}.snv_union.atomized.vcf.gz",
        Q_VCF_multiallelic = "fillout/input/{patient_name}.snv_union.multiallelic.vcf.gz",
        REF_GENOME = config['reference_info']['reference_genome'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    log: 
        "fillout/logs/{patient_name}-patient-wide_genotype.log"
    shell:
        """
        TIMESTAMP=[`date`]
        echo $TIMESTAMP >> {log} 
        echo 'generating candidate alleles from filtered, combined, normed Mutect2 VCF.' >> {log} 

        bgzip -c {input.patient_union_vcf} > {input.patient_union_vcf}.gz && \
        tabix {input.patient_union_vcf}.gz

        # @HZ 03/21/2023 atomize MNPs
        echo 'atomizing MNPs...' >> {log} && \
        {params.BCFTOOLS_EXEC} norm \
            -Oz \
            -o {params.Q_VCF_atomized} \
            --atomize \
            {input.patient_union_vcf}.gz >> {log}

        # @HZ 09/18/2022 get a condensed form for Python use later
        {params.BCFTOOLS_EXEC} query \
            -f '%CHROM:%POS:%REF/%ALT\n' \
            {params.Q_VCF_atomized} > \
            {output.CANDIDATE_ALLELE_for_py}

        # @HZ 08/28/2022: merge multiallelic sites 
        echo 'merging multiallelic sites...' >> {log} && \
        {params.BCFTOOLS_EXEC} norm \
            -Oz \
            -o {params.Q_VCF_multiallelic} \
            --multiallelics + \
            {input.patient_union_vcf}.gz >> {log} 

        {params.BCFTOOLS_EXEC} query \
            -f'%CHROM\t%POS\t%REF,%ALT\n' \
            {params.Q_VCF_multiallelic} | \
        bgzip -c > {output.CANDIDATE_ALLELE} && \
        tabix -s1 -b2 -e2 {output.CANDIDATE_ALLELE} && \
        echo 'finished writing candidate alleles file.' >> {log}

        # cleanup 
        rm {input.patient_union_vcf}.gz {input.patient_union_vcf}.gz.tbi
        rm {params.Q_VCF_atomized} {params.Q_VCF_multiallelic}
        
        echo 'next, bcftools mpileup with:' >> {log} && \
        echo '- Candidate alleles FILE: {output.CANDIDATE_ALLELE}' >> {log} && \
        echo '- GENOME FILE: {params.REF_GENOME}' >> {log} 
        """