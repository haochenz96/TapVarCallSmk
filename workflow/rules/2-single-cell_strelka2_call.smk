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