import pandas as pd
import argparse
import datetime
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--metadata', required=True)
parser.add_argument('--multifasta')
parser.add_argument('--out-dir', help='Directory to save swell files', required=True)
parser.add_argument('-f', '--multifasta-graph')
parser.add_argument('-s', '--start-date')
parser.add_argument('-e', '--end-date')

print("Reading metadata into DataFrame...", end=" ", flush=True)
metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
print("done.")

if args.start_date or args.end_date:
    print("Filtering metadata by given start/end dates...", end=" ", flush=True)
    metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
    metadata = metadata.loc[mask]
    print("done.")

prefix = f"swell_{args.start_date}_{args.end_date}"
filtered_metadata_path = f"{args.out_dir}/{prefix}_metadata.tsv"
print("Saving filtered metadata to tsv in out-dir...", end=" ", flush=True)
metadata.to_csv(filtered_metadata_path, index=False, sep='\t')
print("done.")

if args.multifasta:
    print("Running swell on multifasta...", end=" ", flush=True)
    swell_multifasta = subprocess.run(['swell', 'separate-fasta', args.multifasta], capture_output=True)
    print("done.")
    if args.multifasta_graph:
        print("Generating fasta graphs...", end=" ", flush=True)
        subprocess.run(['Rscript', args.multifasta_graph, swell_multifasta, filtered_metadata_path])
        print("done.")

print("Removing temporary files...", end=" ", flush=True)
subprocess.run(['rm', filtered_metadata_path])
print("done.")