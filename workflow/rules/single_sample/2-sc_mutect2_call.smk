rule step2_sc_mutect2_call_and_add_filters:
    # scattered by single cell
    input:
        sc_bam = "{sample_name}/1-sc_bams/{sample_name}_{cell_num_index}.bam",
    output:  
        vcf = "{sample_name}/2-sc_mutect2_call/m2_sc_vcfs/{sample_name}_{cell_num_index}_somatic_m2.vcf.gz",
        stats = "{sample_name}/2-sc_mutect2_call/m2_sc_vcfs/{sample_name}_{cell_num_index}_somatic_m2.vcf.gz.stats",
        vcf_filter_added = "{sample_name}/2-sc_mutect2_call/m2_sc_vcfs_filter_added/{sample_name}_{cell_num_index}_somatic_m2_filter_added.vcf.gz",
    params:
        REF = config['reference_info']['reference_genome'],
        PANEL_INSERT_FILE = config['reference_info']['panel_insert_file'],
        #GR = config['mutect2']['germline_resource'], # use germline resource as prior prob that normal sample carries an allele
        mrpas = config['mutect2']['mrpas'],
        FILTERM2_OPS = parse_for_shell_args(config['mutect2']['filterm2_ops'])
    retries: 3
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = 8000,
        time_min = lambda wildcards, attempt: attempt * 59 + (attempt-1) * 179, # <------ rate-limiting factor for some cases
    conda:
        "../../envs/mutect2.yaml"
    shell:
    	"sc_bai={input.sc_bam}.bai; "
    	"if ! [ -f '$sc_bai' ]; then samtools index {input.sc_bam}; fi; "
        "gatk Mutect2 "
        "--java-options '-XX:-CreateCoredumpOnCrash' "
        "--max-reads-per-alignment-start {params.mrpas} "
        "-R {params.REF} "
        "-L {params.PANEL_INSERT_FILE} "
        "-I {input.sc_bam} "
        "-O {output.vcf}; "
        "gatk FilterMutectCalls "
        "-R {params.REF} "
        "-V {output.vcf} "
        "-O {output.vcf_filter_added} "
        "{params.FILTERM2_OPS} "
