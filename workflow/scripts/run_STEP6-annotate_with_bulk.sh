# @HZ
# This scripts subsets a VCF file with a given bed file

#!/bin/bash 
#BSUB -J bulk-annotate
#BSUB -sla CMOMEM -R highmem                                                  
#BSUB -n 4                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=4]    # expected resorce consumption for memory
#BSUB -W 2:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/OUTPUTS_from_m2_f/bulk-annotate.stdout
#BSUB -eo /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/OUTPUTS_from_m2_f/bulk-annotate.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate mosaic-custom

python '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/scripts/STEP6-annotate_h5_with_bulk.py' \
    --sample_name RA17_37-39_1 \
    --input_h5 /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/OUTPUTS_from_m2_f/RA17_37-39_1_DNA_CNV_m2_f.h5 \
    --bulk_info_yaml /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/reference/RA17_37-39_1_matched_bulk_info.yaml \
    --output_dir /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/OUTPUTS_from_m2_f \
    --log_file /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA17_37-39_1/OUTPUTS_from_m2_f/RA17_37-39_1.bulk_annot.log