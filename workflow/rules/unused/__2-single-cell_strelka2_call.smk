rule sc_strelk2_configure:
    # (1) create Strelka2-SC configuration script
    # note we are using Strelka2 in germline mode
    # scattered by single cell
    input:
        SC_BAM = "{sample_name}/sc_bams/{sample_name}_{cell_num_index}.bam",
        REF = config['reference_info']['reference_genome'],
        AMPLICON_FILE = config['reference_info']['panel_amplicon_file_gz'],
    output:  
        RUN_SCRIPT = "{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/runWorkflow.py",
    params:
        RUNDIR = "{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/"
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    conda:
        "../envs/strelka2.yaml"
    shell:
    	"configureStrelkaGermlineWorkflow.py "
        "--bam {input.SC_BAM} "
        "--referenceFasta {input.REF} "
        "--runDir {params.RUNDIR} "

rule sc_strelka2_exec:
    # (2) execute Strelka2-SC configuration script
    # note we are using Strelka2 in germline mode
    # scattered by single cell
    input:
        RUN_SCRIPT = "{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/runWorkflow.py",
    output:  
        VCF = "{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/results/variants/variants.vcf.gz",
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        "../envs/strelka2.yaml"
    shell:
        "{input.RUN_SCRIPT} -m local -j {threads}"

# rule combined_strelka2_exec:
#     # (1) create combined_strelka2 configuration script
#     # note we are using Strelka2 in germline mode
#     input:
#         SC_BAMs = get_step1_sc_bams_by_sample,
#         REF = config['reference_info']['reference_genome'],
#         AMPLICON_FILE = config['reference_info']['panel_amplicon_file_gz'],
#     output:  
#         RUN_SCRIPT = "{sample_name}/sc_strelka2_call/runWorkflow.py",
#     params:
#         RUNDIR = "{sample_name}/sc_strelka2_call/"
#     threads: lambda wildcards, attempt: attempt * 2,
#     resources: 
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time_min = lambda wildcards, attempt: attempt * 29,
#     conda:
#         "../envs/strelka2.yaml"
#     shell:
#     	"configureStrelkaGermlineWorkflow.py "
#         "--bam {input.SC_BAM} "
#         "--referenceFasta {input.REF} "
#         "--callRegions {input.AMPLICON_FILE} "
#         "--runDir {params.RUNDIR} "

# rule combined_strelka2_exec:
#     # (2) execute combined_strelka2 configuration script
#     # note we are using Strelka2 in germline mode
#     input:
#         RUN_SCRIPT = "{sample_name}/sc_strelka2_call/runWorkflow.py",
#     output:  
#         VCF = "{sample_name}/sc_strelka2_call/results/variants/variants.vcf.gz",
#     threads: lambda wildcards, attempt: attempt * 12,
#     resources: 
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time_min = lambda wildcards, attempt: attempt * 599,
#     conda:
#         "../envs/strelka2.yaml"
#     shell:
#         "{input.RUN_SCRIPT} -m local -j {threads}"