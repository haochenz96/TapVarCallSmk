#!/bin/bash

# @HZ 02/10/2022
# split a single cells.bam output by the Tapestri pipeline

# ----- split single cell bam -----
module load samtools

while getopts i:o:n: flag
do
    case "${flag}" in
        i) BAM=${OPTARG};;
        o) sc_BAM_dir=${OPTARG};;
        n) sample_name=${OPTARG};;
    esac
done
echo "[INFO] cells.bam --- ${BAM}";
echo "[INFO] sc_bam_dir ---  ${sc_BAM_dir}";
echo "[INFO] sample_name --- ${sample_name}";

mkdir -p $sc_BAM_dir

num_cells=$(samtools view -H $BAM | grep '^@RG' | wc -l)
if [ $(ls $sc_BAM_dir/*.bam | wc -l) != $num_cells ]
then
    echo "[INFO] Single cell BAMs do not exist. Splitting single cell BAMs"
    cd $sc_BAM_dir # generate all the single cell bames here
    # -f flag defines the output filename format
    samtools split $BAM -f "${sample_name}_%!.%." # the string format will make it like cell_barcode.bam
    # index

    echo "[INFO] Splitting single cell BAM done. Iterate and index next."
    for tumor_bam in $(ls ./${sample_name}*.bam)
    do
        samtools index $tumor_bam
    done

else
    echo "[INFO] Single cell BAMs already exist! Making index files if necessary"
    for tumor_bam in $(ls $sc_BAM_dir/*.bam):
    do
        if ! [ -f "${tumor_bam/.bam/.bai}" ]; then samtools index $tumor_bam; fi
    done
fi

