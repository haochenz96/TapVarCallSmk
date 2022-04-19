# (1) get single-cell barcodes for each sample
sample_barcode_map = [] # <--------------------------------------------- global
for sample_i in sample_names:
    bars = [] # <------------------------------------------------------- global
    part_1_output = working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    with pysam.AlignmentFile(part_1_output, "rb") as cells_bam:
        for i in cells_bam.header['RG']:
            #bars.append(i['SM'])
    sample_barcode_map[sample_i] = bars

rule split_sc_bams:
    # scatter per sample, per cell
    input:
        CELLS_BAM = '{sample_name}/tap_pipeline_output/results/bam/{sample_name}.tube1.cells.bam',
        #barcodes = sample_barcode_map[wildcards.sample_name]
    output:
        SC_BAM = "{sample_name}/sc_bams/{sample_name}_{cell_barcode}.bam",
        # sc_bai = expand('sc_bams/{sample_name}_{cell_barcode}.bai', cell_barcode = bars, sample_name = sample_name),
    conda:
        "../envs/samtools.yaml"
    threads: lambda wildcards, attempt: 2**(attempt-1),
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        # -- trial 1 --
        "samtools view -b -r {wildcards.cell_barcode} {input.CELLS_BAM} > {output.SC_BAM} && "
        "samtools index {output.SC_BAM} "

