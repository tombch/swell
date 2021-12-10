import pandas as pd
import argparse
import pathlib
import shutil
import datetime
import subprocess

# Parts of this program were written by Sam Wilkinson and have been modified for usage with swell

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--latest-dir")
parser.add_argument("-o", "--out-dir", help="Directory to save fasta/metadata files", type=pathlib.Path)
parser.add_argument("-g", "--graph-path")
parser.add_argument("--swell-bam", help="", action="store_true")
parser.add_argument("-b", "--bed")
parser.add_argument("-s", "--start-date")
parser.add_argument("-e", "--end-date")
args = parser.parse_args()

print("Reading metadata into DataFrame... ", end="")
metadata_path = f"{args.latest_dir}majora.metadata.matched.tsv"
metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
print("done.")

prefix = f"swell_data_{args.start_date}_{args.end_date}"
out_fasta_path = f"{args.out_dir}/{prefix}.fasta"
out_bam_data_path = f"{args.out_dir}/{prefix}_bam_data.tsv"
out_metadata_path = f"{args.out_dir}/{prefix}_metadata.tsv"

metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
if args.start_date or args.end_date:
    print("Filtering metadata by given start/end dates... ", end="")
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
    metadata = metadata.loc[mask]
    print("done.")
print("Saving metadata to tsv in out-dir... ", end="")
metadata.to_csv(out_metadata_path, index=False, sep='\t')
print("done.")

print("Writing fasta data to multifasta (and optional swell of bam data to a tsv)... ", end="")
with open(out_fasta_path, "w") as out_fasta, open(out_bam_data_path, "w") as out_bam:
    for i, (index, row) in enumerate(metadata.iterrows()):
        fasta_path = f"{args.latest_dir}fasta/{row['central_sample_id']}.{row['run_name']}.climb.fasta"
        with open(fasta_path) as fasta:
            for line in fasta:
                out_fasta.write(line)
        if args.swell_bam:
            bam_path = f"{args.latest_dir}alignment/{row['central_sample_id']}.{row['run_name']}.climb.bam"
            # index_path = f"{args.latest_dir}alignment/{row['central_sample_id']}.{row['run_name']}.climb.bam.bai"
            # bam_name = f"{prefix}_{row['central_sample_id']}.bam"
            # index_name = f"{prefix}_{row['central_sample_id']}.bam.bai"
            # shutil.copyfile(bam_path, f"{args.outdir}/{bam_name}")
            # shutil.copyfile(index_path, f"{args.outdir}/{index_name}")
            swell_bam_out = subprocess.run(['swell', 'bam', bam_path, '--ref', 'NC_045512', 'NC045512', 'MN908947.3', '--bed', args.bed], capture_output=True)
            header, data = (swell_bam_out.stdout.decode('utf8')[:-1]).split('\n')
            if i == 0:
                out_bam.write(f"{header}\n")
                out_bam.write(f"{data}\n")
            else:
                out_bam.write(f"{data}\n")
print("done.")
print("Running swell on multifasta... ", end="")
swell_fasta_out = subprocess.run(['swell', 'separate-fasta', out_fasta_path], capture_output=True)
print("done.")
print("Generating fasta graphs... ", end="")
subprocess.run(['Rscript', args.graph_path, out_metadata_path], input=swell_fasta_out.stdout)
print("done.")
subprocess.run(['rm', out_fasta_path, out_metadata_path])