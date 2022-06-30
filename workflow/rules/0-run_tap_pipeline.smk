# rule rename_barcode_to_num_idx:
#     # scatter per sample
#     input:
#         CELLS_BAM = expand('{sample_name}/tap_pipeline_output/results/bam/{sample_name}.tube1.cells.bam', sample_name = sample_names),
#     output:
#         rename_confirm = '{sample_name}/rename_confirm.txt'
#     params:
#         rename_script = '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/scripts/rename_file_barcode_to_num_idx.sh',
#         bars_map = lambda wildcards: sample_barcode_maps[wildcards.sample_name],
#     run:
#         dir_to_search = f'{wildcards.sample_name}'
#         bar_to_num_map = dict((v,k) for k,v in bars_map.items()) # swap key and value in bars_map
#         for barcode in bar_to_num_map:
#             print('processing barcode: ' + barcode)
#             num_idx = bar_to_num_map[barcode]
#             subprocess.run(['bash', params.rename_script, dir_to_search, barcode, num_idx], capture_output=True)
#         if True:
#             os.system(f'touch {output.rename_confirm}')
        