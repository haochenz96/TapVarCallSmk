rule mutect2_sc_f_pass2:
    # scattered by single cell
    # this step, even at single-cell level, requires lots of resource. Might need some sample-specific hyperparameter tuning.
    input:
        sc_bam = "{sample_name}/sc_bams/{sample_name}_{cell_barcode}.bam",
        q_vcf = "{sample_name}/bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output: 
    	#vcf="$outpath/mutect2/__m2_sc_vcfs/${sample_name}_${cell_name}_somatic_m2.vcf.gz"
        sc_vcf = "{sample_name}/mutect2_sc_f_pass2/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
        stats = "{sample_name}/mutect2_sc_f_pass2/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz.stats",
        sc_vcf_filter_added = "{sample_name}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz",
    params:
        REF = config['reference_info']['reference_genome'],
        mrpas = config['mutect2']['mrpas'],
        FILTERM2_OPS = parse_for_shell_args(config['mutect2']['filterm2_ops']),
    threads: lambda wildcards, attempt: attempt * 2
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 359,
    conda:
        "../envs/var-calling.yaml"
    shell:
    	"sc_bai='{input.sc_bam}.bai'; "
    	"if ! [ -f '$sc_bai' ]; then samtools index {input.sc_bam}; fi; "
        "gatk Mutect2 "
        "--java-options '-XX:-CreateCoredumpOnCrash' "
        "-alleles {input.q_vcf} "
        "-L {input.q_vcf} "
        "--genotype-filtered-alleles "
        "--max-reads-per-alignment-start {params.mrpas} "
        "-R {params.REF} "
        "-I {input.sc_bam} "
        "-O {output.sc_vcf}; "
        "gatk FilterMutectCalls "
        "-R {params.REF} "
        "-V {output.sc_vcf} "
        "-O {output.sc_vcf_filter_added} "
        "{params.FILTERM2_OPS} "