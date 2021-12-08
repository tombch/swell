import pandas as pd
import argparse
import pathlib
import shutil
import datetime

# This program was written by Sam Wilkinson and has been modified for usage before swell

parser = argparse.ArgumentParser()
parser.add_argument("--metadata-path")
parser.add_argument("--latest-dir")
parser.add_argument("-o", "--outdir", help="Directory to save fasta/metadata files", type=pathlib.Path, metavar="")
parser.add_argument("-p", "--prefix", help="Prefix for outfiles", type=str, metavar="")
parser.add_argument("-b", "--get_bam", help="Grab BAM files alongside consensus/metadata", action="store_true")
parser.add_argument("-m", "--metadata", help="Export metadata as TSV file", action="store_true")
parser.add_argument("-s", "--start-date")
parser.add_argument("-e", "--end-date")
args = parser.parse_args()

metadata = pd.read_csv(args.metadata_path, sep="\t", low_memory=False)
metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
if args.start_date or args.end_date:
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
    metadata = metadata.loc[mask]

out_fasta_path = f"{args.outdir}/{args.prefix}.fasta"
with open(out_fasta_path, "w") as out_fasta:
    for index, row in metadata.iterrows():
        fasta_path = f"{args.latest_dir}fasta/{row['central_sample_id']}.{row['run_name']}.climb.fasta"
        with open(fasta_path) as fasta:
            for line in fasta:
                out_fasta.write(line)

    # if args.get_bam:
    #     bam_name = f'{args.prefix}_{row["central_sample_id"]}.bam'
    #     index_name = f'{args.prefix}_{row["central_sample_id"]}.bam.bai'
    #     shutil.copyfile(
    #         "%salignment/%s.%s.climb.bam"
    #         % (args.latest_dir, row["central_sample_id"], row["run_name"]),
    #         f"{args.outdir}/{bam_name}",
    #     )
    #     shutil.copyfile(
    #         "%salignment/%s.%s.climb.bam.bai"
    #         % (args.latest_dir, row["central_sample_id"], row["run_name"]),
    #         f"{args.outdir}/{index_name}",
    #     )

if args.metadata:
    metadata.to_csv(f"{args.outdir}/{args.prefix}.tsv", index=False, sep='\t')