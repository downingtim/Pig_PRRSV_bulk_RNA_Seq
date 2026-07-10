#!/bin/bash

# Directory containing the FASTP files  ->  srun --pty --export=ALL --mem=500G $@ /bin/bash 
fastp_files_dir="fastp/*/*/*"

base_command="kraken2"
db_path="/home/share/"
threads=3
output_dir_base="KRAKEN_FILES"
output_dir_base2="KRAKEN_REPORT"

# Get a list of sample names by listing the files in the directory and extracting the base names
sample_names=$(ls $fastp_files_dir | grep -E '_[12].fq.gz$' | sed -E 's/_([12]).fq.gz$//' | sort | uniq)

# Iterate over sample names and construct the command for each
for sample_name in $sample_names; do
    sample2="${sample_name##*/}"    # Construct output directories
    output_dir="${output_dir_base}/${sample2}"
    output_dir2="${output_dir_base2}/${sample2}"
    command="${base_command} --db ${db_path} --threads ${threads} --report ${output_dir2} --use-names --output ${output_dir} --paired ${sample_name}_1.fq.gz ${sample_name}_2.fq.gz " 
#   command="${base_command} --db ${db_path} --threads ${threads} --report ${output_dir2} --use-names --output ${output_dir} --paired ${sample_name}_1.fq.gz ${sample_name}_2.fq.gz " 
    echo "${command}"
done
