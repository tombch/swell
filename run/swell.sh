#!/usr/bin/bash
out_dir=$1
latest_dir=$2
start_date=$3
end_date=$4
prefix="swell_${start_date}_${end_date}"
metadata_path="${latest_dir}majora.metadata.matched.tsv"
python3 get_seqs_metadata.py --outdir $out_dir --prefix $prefix --metadata-path $metadata_path --latest-dir $latest_dir -s $start_date -e $end_date --metadata
out_fasta_path="${out_dir}/${prefix}.fasta"
out_metadata_path="${out_dir}/${prefix}.tsv"
swell separate-fasta $out_fasta_path | Rscript ../separate_fasta_plot.R $out_metadata_path