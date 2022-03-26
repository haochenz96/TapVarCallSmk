rule split_sc_bams:
    input:
        CELLS_BAM = part_1_output,
    output:
        #sc_bams_dir = "sc_bams",
        sc_bam = expand('sc_bams/{sample_name}_{cell_barcode}.bam', cell_barcode = bars, sample_name = sample_name),
        # sc_bai = expand('sc_bams/{sample_name}_{cell_barcode}.bai', cell_barcode = bars, sample_name = sample_name),
    params:
        sc_bams_dir = "sc_bams",
        sample_name = sample_name,
        script = scripts_dir / 'split_cellbam.sh'
    conda:
        "../envs/samtools.yaml"
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = lambda wildcards, attempt: attempt * 59,
    shell:
        "mkdir -p {params.sc_bams_dir}; "
        "cd {params.sc_bams_dir}; "
        "samtools split {input.CELLS_BAM} -f '{params.sample_name}_%!.%.'"
        # "mkdir -p {params.sc_bams_dir}; "
        # "bash {params.script} "
        # "-i {input.CELLS_BAM} "
        # "-o {params.sc_bams_dir} "
        # "-n {params.sample_name} "


# rule get_sc_bam_mapping:
#     input:
#         sc_bam = expand('sc_bams/{sample_name}_{cell_barcode}.bam', cell_barcode = bars, sample_name = sample_name),
#         # sc_bai = expand('sc_bams/{sample_name}_{cell_barcode}.bai', cell_barcode = bars, sample_name = sample_name),
#     output:
#         barcode_map = working_dir / 'sc_barcode_map.txt'
#     params:
#         sc_bam_dir = 'sc_bams'
#     run:
#         get_single_cell_name_mapping(
#             params.sc_bam_dir, 
#             output.barcode_map
#             )

rule mutect2_sc_pass1:
    # scattered by single cell
    input:
        sc_bam = "sc_bams/{sample_name}_{cell_barcode}.bam",
    output: 
    	#vcf="$outpath/mutect2/__m2_sc_vcfs/${sample_name}_${cell_name}_somatic_m2.vcf.gz"
        vcf = "mutect2_sc_pass1/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
        stats = "mutect2_sc_pass1/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz.stats",
        vcf_filter_added = "mutect2_sc_pass1/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz",
    params:
        REF = config['reference_info']['reference_genome'],
        AMPLICON_FILE = config['reference_info']['panel_amplicon_file'],
        GR = config['mutect2']['germline_resource'], # use germline resource as prior prob that normal sample carries an allele
        mrpas = config['mutect2']['mrpas'],
        FILTERM2_OPS = parse_for_shell_args(config['mutect2']['filterm2_ops'])
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/var-calling.yaml"
    shell:
    	"sc_bai={input.sc_bam}.bai; "
    	"if ! [ -f '$sc_bai' ]; then samtools index {input.sc_bam}; fi; "
        "gatk Mutect2 "
        "--max-reads-per-alignment-start {params.mrpas} "
        "--germline-resource {params.GR} "
        "-R {params.REF} "
        "-L {params.AMPLICON_FILE} "
        "-I {input.sc_bam} "
        "-O {output.vcf}; "
        "gatk FilterMutectCalls "
        "-R {params.REF} "
        "-V {output.vcf} "
        "-O {output.vcf_filter_added} "
        "{params.FILTERM2_OPS} "

# rule FilterMutectCall:
# # scattered by single cell
#     input:
#         vcf = "mutect2_sc_pass1/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
#     output: 
#         vcf_filter_added = "mutect2_sc_pass1/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz",
#     params:
#         REF = config['reference_info']['reference_genome'],
#         # AMPLICON_FILE = config['reference_info']['panel_amplicon_file'],
#         FILTERM2_OPS = parse_for_shell_args(config['mutect2']['filterm2_ops'])
#     resources: 
#         mem_mb = lambda wildcards, attempt: attempt * 2000,
#         time_min = 10
#     conda:
#         "../envs/var-calling.yaml"
#     shell:
#         "gatk FilterMutectCalls "
#         "-R {params.REF} "
#         "-V {input.vcf} "
#         "-O {output.vcf_filter_added} "
#         "{params.FILTERM2_OPS}"


