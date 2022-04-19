# (1) get single-cell barcodes for each sample
sample_barcode_maps = {} # <--------------------------------------------- global
for sample_i in sample_names:
    bars_map = {} # <------------------------------------------------------- global
    part_1_output = working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    with pysam.AlignmentFile(part_1_output, "rb") as cells_bam:
        # @HZ 04/19/2022 create numerical index for each barcode
        cell_count = 0
        for i in cells_bam.header['RG']:
            cell_count += 1
            bars_map[f'cell_{cell_count}'] = i['SM'] # <-------------------- numerical index to cell barcode map
            #bars.append(i['SM'])

    sample_barcode_maps[sample_i] = bars

rule write_barcode_map:
    # scatter per sample
    input:
        CELLS_BAM = expand('{sample_name}/tap_pipeline_output/results/bam/{sample_name}.tube1.cells.bam', sample_name = sample_names),
    output:
        BARCODE_MAP = expand('{sample_name}/references/{sample_name}.barcode_map.txt', sample_name = sample_names),
    threads: lambda wildcards, attempt: 2**(attempt-1),
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 29,
    run:
        # for each sample, write a barcode map
        for sample_i in sample_names:
            bars_map = sample_barcode_maps[sample_i]
            bars_map_df = pd.DataFrame.from_dict(bars_map, orient='index', columns = ['cell_barcode'])
            bars_map_df.index.name = 'num_index'
            bars_map_df.to_csv(f'{sample_i}/references/{sample_i}.barcode_map.txt', sep='\t', index=True)

rule split_sc_bams:
    # scatter per sample, per cell
    input:
        CELLS_BAM = '{sample_name}/tap_pipeline_output/results/bam/{sample_name}.tube1.cells.bam',
        #barcodes = sample_barcode_map[wildcards.sample_name]
    output:
        SC_BAM = "{sample_name}/sc_bams/{sample_name}_{cell_num_index}.bam",
        # sc_bai = expand('sc_bams/{sample_name}_{cell_barcode}.bai', cell_barcode = bars, sample_name = sample_name),
    params:
        cell_barcode = sample_barcode_maps[{wildcards.sample_name}][{wildcards.cell_num_index}],
    conda:
        "../envs/samtools.yaml"
    threads: lambda wildcards, attempt: 2**(attempt-1),
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        # -- trial 1 --
        "samtools view -b -r {params.cell_barcode} {input.CELLS_BAM} > {output.SC_BAM} && "
        "samtools index {output.SC_BAM} "