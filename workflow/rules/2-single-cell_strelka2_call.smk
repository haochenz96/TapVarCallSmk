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

rule strelk2_sc_configure:
    # scattered by single cell
    input:
        SC_BAM = "{sample_name}/sc_bams/{sample_name}_{cell_barcode}.bam",
    output:  
        RUN_SCRIPT = "{sample_name}/strelka2_sc/{sample_name}_{cell_barcode}_germline_s2.runWorkflow.py",
    params:
        REF = config['reference_info']['reference_genome'],
        AMPLICON_FILE = config['reference_info']['panel_amplicon_file'],
        RUNDIR = "strelka/{patient}/{sample}/"
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/var-calling.yaml"
    shell:
    	"configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal} "
        "--tumorBam {input.tumor} "
        "--referenceFasta {input.ref} "
        "--runDir {params.rundir} "

rule strelka2_sc_exec:
    # run each strelka configureation script
    input:
        runpy = "strelka/{patient}/{sample}/runWorkflow.py"
    output:  
        VCF = "{sample_name}/strelka2_sc/s2_sc_vcfs/{sample_name}_{cell_barcode}_germline_s2.vcf.gz",
    params:
        cores = config['strelka2']['cores']
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "{input.runpy} -m local -j {params.cores}"