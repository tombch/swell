#!/usr/bin/bash
metadata_path=/cephfs/covid/bham/artifacts/published/latest/majora.metadata.matched.tsv
start_date=$1
end_date=$2
python3 get_seqs_metadata.py -s $start_date -e $end_date | swell separate-fasta - | Rscript separate_fasta_plot.R $metadata_path