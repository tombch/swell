import pandas as pd
import argparse
import pathlib
import shutil
import datetime

# Code written by Sam Wilkinson

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outdir", help="Directory to save fasta/metadata files", type=pathlib.Path, metavar="")
parser.add_argument("-p", "--prefix", help="Prefix for outfiles", type=str, metavar="")
parser.add_argument("-b", "--get_bam", help="Grab BAM files alongside consensus/metadata", action="store_true")
parser.add_argument("-m", "--metadata", help="Export metadata as TSV file", action="store_true")
parser.add_argument("-s", "--start-date")
parser.add_argument("-e", "--end-date")
args = parser.parse_args()

metadata = pd.read_csv("/cephfs/covid/bham/artifacts/published/latest/majora.metadata.matched.tsv", sep="\t", low_memory=False)
metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
if args.start_date or args.end_date:
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
    metadata = metadata.loc[mask]

latest_dir = "/cephfs/covid/bham/artifacts/published/latest/"

if args.outdir and args.prefix:
    out_fasta = open("%s/%s.fasta" % (args.outdir, args.prefix), "w")

for index, row in metadata.iterrows():
    fasta_path = "%sfasta/%s.%s.climb.fasta" % (
        latest_dir,
        row["central_sample_id"],
        row["run_name"],
    )
    with open(fasta_path) as f:
        if args.outdir and args.prefix:
            for line in f:
                out_fasta.write(line)
        else:
            for line in f:
                # Use for piping only!
                print(line)

    # if args.get_bam:
    #     bam_name = f'{args.prefix}_{row["central_sample_id"]}.bam'
    #     index_name = f'{args.prefix}_{row["central_sample_id"]}.bam.bai'
    #     shutil.copyfile(
    #         "%salignment/%s.%s.climb.bam"
    #         % (latest_dir, row["central_sample_id"], row["run_name"]),
    #         f"{args.outdir}/{bam_name}",
    #     )
    #     shutil.copyfile(
    #         "%salignment/%s.%s.climb.bam.bai"
    #         % (latest_dir, row["central_sample_id"], row["run_name"]),
    #         f"{args.outdir}/{index_name}",
    #     )

if args.metadata:
    metadata.to_csv("%s/%s.csv" % (args.outdir, args.prefix), index=False)