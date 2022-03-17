rule mutect2_sc_f_pass2:
    # scattered by single cell
    input:
        sc_bam = "sc_bams/{sample_name}_{cell_barcode}.bam",
        q_vcf = "bcf_pass1/combined_vcf/{sample_name}-combined_VCF_filter_v3.prev_filtered.vcf.gz",
    output: 
    	#vcf="$outpath/mutect2/__m2_sc_vcfs/${sample_name}_${cell_name}_somatic_m2.vcf.gz"
        sc_vcf = "mutect2_sc_f_pass2/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
    params:
        REF = config['reference_info']['reference_genome'],
        mrpas = config['mutect2']['mrpas']
    threads: lambda wildcards, attempt: attempt * 2
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = 179
    conda:
        "../envs/var-calling.yaml"
    shell:
    	"sc_bai='{input.sc_bam}.bai'; "
    	"if ! [ -f '$sc_bai' ]; then samtools index {input.sc_bam}; fi; "
        "gatk Mutect2 "
        "-alleles {input.q_vcf} "
        "-L {input.q_vcf} "
        "--genotype-filtered-alleles "
        "--max-reads-per-alignment-start {params.mrpas} "
        "-R {params.REF} "
        "-I {input.sc_bam} "
        "-O {output.sc_vcf} "