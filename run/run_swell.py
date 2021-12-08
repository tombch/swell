import pandas as pd
import argparse
import pathlib
import shutil
import datetime
import subprocess

# This program was written by Sam Wilkinson and has been modified for usage with swell

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data-dir")
parser.add_argument("-o", "--out-dir", help="Directory to save fasta/metadata files", type=pathlib.Path, metavar="")
# parser.add_argument("-b", "--get_bam", help="Grab BAM files alongside consensus/metadata", action="store_true")
parser.add_argument("-s", "--start-date")
parser.add_argument("-e", "--end-date")
args = parser.parse_args()

metadata_path = f"{args.data_dir}majora.metadata.matched.tsv"
metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")

prefix = f"swell_data_{args.start_date}_{args.end_date}"
out_fasta_path = f"{args.out_dir}/{prefix}.fasta"
out_metadata_path = f"{args.out_dir}/{prefix}.tsv"
metadata.to_csv(out_metadata_path, index=False, sep='\t')

if args.start_date or args.end_date:
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
    metadata = metadata.loc[mask]

with open(out_fasta_path, "w") as out_fasta:
    for index, row in metadata.iterrows():
        fasta_path = f"{args.data_dir}fasta/{row['central_sample_id']}.{row['run_name']}.climb.fasta"
        with open(fasta_path) as fasta:
            for line in fasta:
                out_fasta.write(line)

        # if args.get_bam:
        #     bam_path = f"{args.data_dir}alignment/{row['central_sample_id']}.{row['run_name']}.climb.bam"
        #     index_path = f"{args.data_dir}alignment/{row['central_sample_id']}.{row['run_name']}.climb.bam.bai"
        #     bam_name = f"{prefix}_{row['central_sample_id']}.bam"
        #     index_name = f"{prefix}_{row['central_sample_id']}.bam.bai"
        #     shutil.copyfile(bam_path, f"{args.outdir}/{bam_name}")
        #     shutil.copyfile(index_path, f"{args.outdir}/{index_name}")

swell_out = subprocess.run(['swell', 'separate-fasta', out_fasta_path], capture_output=True)
subprocess.run(['Rscript', '../separate_fasta_plot.R', out_metadata_path], input=swell_out.stdout)
subprocess.run(['rm', out_fasta_path, out_metadata_path])