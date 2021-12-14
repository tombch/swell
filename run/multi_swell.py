import pandas as pd
import argparse
import subprocess


# This function was written by Heng Li
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


parser = argparse.ArgumentParser()
parser.add_argument('--metadata', required=True, help='Path to sample metadata')
parser.add_argument('--temp-dir', required=True, help='Directory to save temporary files')
parser.add_argument('--multifasta', help='Path to multifasta file')
parser.add_argument('--multifasta-graph', help='Path to graphs that use swell stats generated from the multifasta')
parser.add_argument('-s', '--start-date', help='Minimum sequencing_submission_date for a sample to be included in swell run')
parser.add_argument('-e', '--end-date', help='Maximum sequencing_submission_date for a sample to be included in swell run')
args = parser.parse_args()


print("Reading metadata...", end=" ", flush=True)
metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
print("done.")


if args.start_date or args.end_date:
    print("Filtering metadata by the provided start/end dates...", end=" ", flush=True)
    date_filtered_header_dict = {}
    start_date = pd.Timestamp(args.start_date)
    end_date = pd.Timestamp(args.end_date)
    metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
    
    if args.start_date and args.end_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
        for i, (index, row) in enumerate(metadata.iterrows()):
            if row['sequencing_submission_date'] >= start_date and row['sequencing_submission_date'] <= end_date:
                date_filtered_header_dict[row['fasta_header']] = True
            else:
                date_filtered_header_dict[row['fasta_header']] = False
    elif args.start_date:
        mask = (metadata['sequencing_submission_date'] >= args.start_date)
        for i, (index, row) in enumerate(metadata.iterrows()):
            if row['sequencing_submission_date'] >= start_date:
                date_filtered_header_dict[row['fasta_header']] = True
            else:
                date_filtered_header_dict[row['fasta_header']] = False
    elif args.end_date:
        mask = (metadata['sequencing_submission_date'] <= args.end_date)
        for i, (index, row) in enumerate(metadata.iterrows()):
            if row['sequencing_submission_date'] <= end_date:
                date_filtered_header_dict[row['fasta_header']] = True
            else:
                date_filtered_header_dict[row['fasta_header']] = False
    
    metadata = metadata.loc[mask]
    print("done.")


prefix = f"swell_{args.start_date}_{args.end_date}"
date_filtered_metadata_path = f"{args.temp_dir}/{prefix}_metadata.tsv"
date_filtered_multifasta_path = f"{args.temp_dir}/{prefix}_multifasta.fasta"
metadata.to_csv(date_filtered_metadata_path, index=False, sep='\t')


if args.multifasta:
    if args.start_date or args.end_date:
        with open(args.multifasta) as multifasta, open(date_filtered_multifasta_path, "w") as date_filtered_multifasta:
            heng_iter = readfq(multifasta)
            for name, seq, qual in heng_iter:
                header = (name.split('|'))[0]
                if date_filtered_header_dict[header]:
                    date_filtered_multifasta.write(f">{header}\n{seq}\n")
        multifasta_path = date_filtered_multifasta_path
    else:
        multifasta_path = args.multifasta
    
    print("Running swell on multifasta...", end=" ", flush=True)
    swell_multifasta = subprocess.run(['swell', 'separate-fasta', multifasta_path], capture_output=True)
    print("done.")
    
    if args.multifasta_graph:
        print("Generating graphs of swell multifasta data...")
        subprocess.run(['Rscript', args.multifasta_graph, date_filtered_metadata_path], input=swell_multifasta.stdout)
        print("done.")


print("Removing temporary files...", end=" ", flush=True)
subprocess.run(['rm', date_filtered_metadata_path, date_filtered_multifasta_path])
print("done.")