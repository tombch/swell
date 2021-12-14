import pandas as pd
import argparse
import pathlib
import shutil
import datetime
import subprocess

# Parts of this program were written by Sam Wilkinson and have been modified for usage with swell

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--latest-dir")
parser.add_argument("-o", "--out-dir", help="Directory to save swell files", type=pathlib.Path)
parser.add_argument("--swell-fasta", action="store_true")
parser.add_argument("--swell-bam", action="store_true")
parser.add_argument("--bed")
parser.add_argument("-f", "--fasta-graph-path")
parser.add_argument("-b", "--bam-graph-path")
parser.add_argument("-s", "--start-date")
parser.add_argument("-e", "--end-date")
args = parser.parse_args()

prefix = f"swell_data_{args.start_date}_{args.end_date}"
metadata_path = f"{args.latest_dir}majora.metadata.matched.tsv"
print("Reading metadata into DataFrame...", end=" ", flush=True)
metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
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

out_metadata_path = f"{args.out_dir}/{prefix}_metadata.tsv"
print("Saving metadata to tsv in out-dir...", end=" ", flush=True)
metadata.to_csv(out_metadata_path, index=False, sep='\t')
print("done.")

if args.swell_fasta:
    print("Running swell on fasta files and writing results into out-dir...", end=" ", flush=True)
    out_fasta_data_path = f"{args.out_dir}/{prefix}_fasta_data.tsv"
    with open(out_fasta_data_path, "w") as out_fasta_data:
        for i, (index, row) in enumerate(metadata.iterrows()):
            fasta_path = f"{args.latest_dir}fasta/{row['central_sample_id']}.{row['run_name']}.climb.fasta"
            swell_fasta_out = subprocess.run(['swell', 'separate-fasta', fasta_path], capture_output=True)
            header, data = (swell_fasta_out.stdout.decode('utf8')[:-1]).split('\n')
            if i == 0:
                out_fasta_data.write(f"{header}\n")
                out_fasta_data.write(f"{data}\n")
            else:
                out_fasta_data.write(f"{data}\n")
    print("done.")

if args.swell_bam:
    print("Running swell on bam files and writing results into out-dir...", end=" ", flush=True)
    out_bam_data_path = f"{args.out_dir}/{prefix}_bam_data.tsv"
    with open(out_bam_data_path, "w") as out_bam_data:
        for i, (index, row) in enumerate(metadata.iterrows()):
            bam_path = f"{args.latest_dir}alignment/{row['central_sample_id']}.{row['run_name']}.climb.bam"
            swell_bam_out = subprocess.run(['swell', 'bam', bam_path, '--ref', 'NC_045512', 'NC045512', 'MN908947.3', '--bed', args.bed], capture_output=True)
            header, data = (swell_bam_out.stdout.decode('utf8')[:-1]).split('\n')
            if i == 0:
                out_bam_data.write(f"{header}\n")
                out_bam_data.write(f"{data}\n")
            else:
                out_bam_data.write(f"{data}\n")
    print("done.")

if args.fasta_graph_path:
    print("Generating fasta graphs...", end=" ", flush=True)
    subprocess.run(['Rscript', args.fasta_graph_path, out_fasta_data_path, out_metadata_path])
    print("done.")

if args.bam_graph_path:
    print("Generating bam graphs...", end=" ", flush=True)
    subprocess.run(['Rscript', args.bam_graph_path, out_bam_data_path, out_metadata_path])
    print("done.")

print("Removing temporary data files...", end=" ", flush=True)
subprocess.run(['rm', out_fasta_data_path, out_bam_data_path, out_metadata_path])
print("done.")