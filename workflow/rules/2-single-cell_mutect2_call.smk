rule split_sc_bams:
    input:
        CELLS_BAM = '{sample_name}/tap_pipeline_output/results/bam/{sample_name}.tube1.cells.bam',
        #barcodes = sample_barcode_map[wildcards.sample_name]
    output:
        #sc_bams_dir = "sc_bams",
        sc_bam = "{sample_name}/sc_bams/{sample_name}_{cell_barcode}.bam",
        # sc_bai = expand('sc_bams/{sample_name}_{cell_barcode}.bai', cell_barcode = bars, sample_name = sample_name),
    # params:
    #     sc_bams_dir = "{sample_name}/sc_bams",
    #     #script = scripts_dir / 'split_cellbam.sh'
    conda:
        "../envs/samtools.yaml"
    threads: lambda wildcards, attempt: 2**(attempt-1),
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        # -- trial 1 --
        "samtools view -b -r {wildcards.cell_barcode} {input.CELLS_BAM} > {output.sc_bam}"
        # -- trial 2 --
        # "mkdir -p {output.sc_bams_dir}; "
        # "cd {output.sc_bams_dir}; "
        # "samtools split {input.CELLS_BAM} -f '{wildcards.sample_name}_%!.%.'"
        # "mkdir -p {params.sc_bams_dir}; "
        # "bash {params.script} "
        # "-i {input.CELLS_BAM} "
        # "-o {params.sc_bams_dir} "
        # "-n {params.sample_name} "

rule mutect2_sc_pass1:
    # scattered by single cell
    input:
        sc_bam = "{sample_name}/sc_bams/{sample_name}_{cell_barcode}.bam",
    output:
         
        vcf = "{sample_name}/mutect2_sc_pass1/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz",
        stats = "{sample_name}/mutect2_sc_pass1/m2_sc_vcfs/{sample_name}_{cell_barcode}_somatic_m2.vcf.gz.stats",
        vcf_filter_added = "{sample_name}/mutect2_sc_pass1/m2_sc_vcfs_filter_added/{sample_name}_{cell_barcode}_somatic_m2_filter_added.vcf.gz",
    params:
        REF = config['reference_info']['reference_genome'],
        AMPLICON_FILE = config['reference_info']['panel_amplicon_file'],
        #GR = config['mutect2']['germline_resource'], # use germline resource as prior prob that normal sample carries an allele
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
        "--java-options '-XX:-CreateCoredumpOnCrash' "
        "--max-reads-per-alignment-start {params.mrpas} "
        "-R {params.REF} "
        "-L {params.AMPLICON_FILE} "
        "-I {input.sc_bam} "
        "-O {output.vcf}; "
        "gatk FilterMutectCalls "
        "-R {params.REF} "
        "-V {output.vcf} "
        "-O {output.vcf_filter_added} "
        "{params.FILTERM2_OPS} "
