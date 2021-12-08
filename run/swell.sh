#!/usr/bin/bash 
latest_dir=$1
start_date=$2
end_date=$3
metadata_path="${latest_dir}majora.metadata.matched.tsv"
python3 get_seqs_metadata.py --metadata-path $metadata_path --latest-dir $latest_dir -s $start_date -e $end_date | swell separate-fasta - | Rscript ../separate_fasta_plot.R $metadata_path