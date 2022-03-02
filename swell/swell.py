import re
import sys
import csv
import pysam
import argparse
import numpy as np
import concurrent.futures
from . import readfq # thanks heng


def clip_tiles(tiles):
    # iterate through tiles and clip
    new_tiles = []
    for tile_index, tile_tuple in enumerate(tiles):
        tile_dict = dict(tile_tuple[2]) # copy the dict god this is gross stuff

        # Clip the start of this window to the end of the last window
        # (if there is a last window)
        if tile_index > 0:
            tile_dict["inside_start"] = tiles[tile_index - 1][2]["end"]

        # Clip the end of this window to the start of the next window
        # (if there is a next window)
        if tile_index < len(tiles) - 1:
            tile_dict["inside_end"] = tiles[tile_index + 1][2]["start"]

        new_tiles.append((tile_tuple[0], tile_tuple[1], tile_dict))
    return new_tiles


def load_scheme(scheme_bed, clip=True):
    tiles_dict = {}
    scheme_fh = open(scheme_bed)
    for line in scheme_fh:
        data = line.strip().split()
        start, end, tile = data[1], data[2], data[3] 
        scheme, tile, side = tile.split("_", 2)

        start = int(start)
        end = int(end)

        if tile not in tiles_dict:
            tiles_dict[tile] = {
                "start": -1,
                "inside_start": -1,
                "inside_end": -1,
                "end": -1,
            }

        if "LEFT" in side.upper():
            if tiles_dict[tile]["start"] == -1:
                tiles_dict[tile]["start"] = start
                tiles_dict[tile]["inside_start"] = end

            if start < tiles_dict[tile]["start"]:
                # push the window region to the leftmost left position
                tiles_dict[tile]["start"] = start
            if end > tiles_dict[tile]["inside_start"]:
                # open the start of the inner window to the rightmost left position
                tiles_dict[tile]["inside_start"] = end

        elif "RIGHT" in side.upper():
            if tiles_dict[tile]["end"] == -1:
                tiles_dict[tile]["end"] = end
                tiles_dict[tile]["inside_end"] = start

            if end > tiles_dict[tile]["end"]:
                # stretch the window out to the rightmost right position
                tiles_dict[tile]["end"] = end
            if start < tiles_dict[tile]["inside_end"]:
                # close the end of the inner window to the leftmost right position
                tiles_dict[tile]["inside_end"] = start
    
    tiles_list = []
    tiles_seen = set([])
    scheme_fh.seek(0)
    for line in scheme_fh:
        data = line.strip().split()
        start, end, tile = data[1], data[2], data[3] 
        scheme, tile, side = tile.split("_", 2)
        tile_tup = (scheme, tile, tiles_dict[tile])
        if tiles_dict[tile]["inside_start"] != -1 and tiles_dict[tile]["inside_end"] != -1 and tile not in tiles_seen:
            tiles_list.append(tile_tup)
            tiles_seen.add(tile)

    tiles_list = sorted(tiles_list, key=lambda x: int(x[1])) # sort by tile number
    if clip: # Default
        new_tiles = clip_tiles(tiles_list)
    else:
        new_tiles = tiles_list

    return new_tiles


def swell_from_fasta(fasta_path):
    if fasta_path == "-":
        fastas = readfq.readfq(sys.stdin)
    else:
        fastas = readfq.readfq(open(fasta_path))
        
    rows = []
    for name, seq, qual in fastas:
        num_seqs = 1
        num_bases = 0
        num_acgt = 0
        num_masked = 0
        num_invalid = 0
        num_ambiguous = 0
        current_gap = 0
        current_ungap = 0
        max_gap = 0
        max_ungap = 0
        prop_acgt = 0
        prop_masked = 0
        prop_invalid = 0
        prop_ambiguous = 0

        for base in seq:
            num_bases += 1

            if base.upper() in 'ACGT':
                num_acgt += 1
                current_ungap += 1
                if current_gap > max_gap:
                    max_gap = current_gap
                current_gap = 0

            elif base.upper() in 'WSMKRYBDHV':
                num_ambiguous += 1
                current_ungap += 1
                if current_gap > max_gap:
                    max_gap = current_gap
                current_gap = 0
            
            elif base.upper() in 'NX':
                num_masked += 1
                current_gap += 1
                if current_ungap > max_ungap:
                    max_ungap = current_ungap
                current_ungap = 0
            
            else:
                num_invalid += 1
                current_gap += 1
                if current_ungap > max_ungap:
                    max_ungap = current_ungap
                current_ungap = 0

        if num_bases > 0:
            prop_acgt = num_acgt / num_bases * 100.0
            prop_masked = num_masked / num_bases * 100.0
            prop_invalid = num_invalid / num_bases * 100.0
            prop_ambiguous = num_ambiguous / num_bases * 100.0
        else:
            prop_invalid = 100.0

        # TODO: is this reliable
        header = (name.split('|'))[0]
        rows.append([fasta_path, header, num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap])

    return ["fasta_path", "header", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], rows


def summarise_swell_from_fasta(fasta_path):
    if fasta_path == "-":
        fastas = readfq.readfq(sys.stdin)
    else:
        fastas = readfq.readfq(open(fasta_path))

    num_seqs = 0
    num_bases = 0
    num_acgt = 0
    num_masked = 0
    num_invalid = 0
    num_ambiguous = 0
    current_gap = 0
    current_ungap = 0
    max_gap = 0
    max_ungap = 0
    prop_acgt = 0
    prop_masked = 0
    prop_invalid = 0
    prop_ambiguous = 0

    rows = []
    for name, seq, qual in fastas:
        num_seqs += 1

        for base in seq:
            num_bases += 1

            if base.upper() in 'ACGT':
                num_acgt += 1
                current_ungap += 1
                if current_gap > max_gap:
                    max_gap = current_gap
                current_gap = 0

            elif base.upper() in 'WSMKRYBDHV':
                num_ambiguous += 1
                current_ungap += 1
                if current_gap > max_gap:
                    max_gap = current_gap
                current_gap = 0
            
            elif base.upper() in 'NX':
                num_masked += 1
                current_gap += 1
                if current_ungap > max_ungap:
                    max_ungap = current_ungap
                current_ungap = 0
            
            else:
                num_invalid += 1
                current_gap += 1
                if current_ungap > max_ungap:
                    max_ungap = current_ungap
                current_ungap = 0

    if num_bases > 0:
        prop_acgt = num_acgt / num_bases * 100.0
        prop_masked = num_masked / num_bases * 100.0
        prop_invalid = num_invalid / num_bases * 100.0
        prop_ambiguous = num_ambiguous / num_bases * 100.0
    else:
        prop_invalid = 100.0

    rows.append([fasta_path, "-", num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap])
    
    return ["fasta_path", "header", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], rows


def swell_from_depth_iter(depth_iterable, depth_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    threshold_counters = {threshold : 0 for threshold in thresholds}
    tile_threshold_counters = {threshold : 0 for threshold in thresholds}
    n_positions = 0
    avg_cov = 0
    n_lines = 0

    tile_starts = [t[2]["inside_start"] for t in tiles] # dont use -1 for 1-pos depth files
    tile_ends = [t[2]["inside_end"] for t in tiles]
    
    tile_data = [[] for t in tiles]
    for line in depth_iterable:
        n_lines += 1
        ref, pos, cov = line.strip().split('\t')
        if sum([g in ref for g in genomes]) != 1:
            continue
        pos = int(pos)
        cov = int(cov)

        # Count positions above threshold
        for threshold in threshold_counters:
            if cov >= threshold:
                threshold_counters[threshold] += 1
        n_positions += 1
        avg_cov = avg_cov + (cov - avg_cov)/n_positions

        # Add coverage value to every open tile at this position
        if tiles:
            for t_i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
                if start <= pos <= end:
                    tile_data[t_i].append(cov)

    tile_vector = []
    for covs in tile_data:
        median_cov = np.median(covs)
        tile_vector.append(median_cov)

        # Count tile means above threshold
        for threshold in threshold_counters:
            if median_cov >= threshold:
                tile_threshold_counters[threshold] += 1

    if min_pos:
        if n_positions < min_pos:
            if min_pos_total_zero and n_lines == 0:
                #TODO This probably only works for depth with -a not -aa?
                sys.stderr.write("[FAIL] BAM has no reads aligned to allowed reference list, but we'll ignore this as seems to have no reads at all.")
            else:
                sys.stderr.write("[FAIL] BAM has fewer than %d positions covered by the allowed reference list.\n" % min_pos)
                sys.exit(3)

    if n_positions > 0:
        threshold_counts_prop = [threshold_counters[x]/n_positions * 100.0 for x in sorted(thresholds)]
    else:
        threshold_counts_prop = [0 for x in sorted(thresholds)]

    if len(tiles) > 0:
        tile_threshold_counts_prop = [tile_threshold_counters[x]/len(tiles) * 100.0 for x in sorted(thresholds)]
    else:
        tile_threshold_counts_prop = [0 for x in sorted(thresholds)]

    if len(tile_vector) == 0:
        tile_vector_str = "-"
    else:
        tile_vector_str = ",".join(["%.2f" % x for x in tile_vector])
    
    rows = []
    rows.append([depth_path.replace(".depth", ""), n_positions, avg_cov] + threshold_counts_prop + tile_threshold_counts_prop + [len(tile_vector), tile_vector_str])
    return ["bam_path", "num_pos", "mean_cov"] + ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)] + ["pc_tiles_medcov_gte%d" % x for x in sorted(thresholds)] + ["tile_n", "tile_vector"], rows 


def swell_from_depth(depth_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    depth_fh = open(depth_path)
    return swell_from_depth_iter(depth_fh, depth_path, tiles, genomes, thresholds, min_pos, min_pos_total_zero)


def swell_from_bam(bam_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    # samtools 1.14
    # pysam 0.18.0 - much better than 0.16
    depth_iterable = (x.group(0)[:-1] for x in re.finditer('.*\n', pysam.depth('-a', bam_path))) # type: ignore
    return swell_from_depth_iter(depth_iterable, bam_path, tiles, genomes, thresholds, min_pos, min_pos_total_zero)


def swell_from_row(args):
    # Assign individual variable names to the arguments
    record, genomes, metadata_headers, thresholds, dp, min_pos, min_pos_total_zero, clip = args

    if record['ref'] and (not record['ref'] in set(genomes)):
        genomes.append(record['ref'])
    tiles = {}
    if record['bed_path']:
        tiles = load_scheme(record['bed_path'], clip)
    
    # Run swell fasta and swell bam on the paths given in the record. Then, return the outputted fields as a tab-separated string
    _, fields = swell_from_fasta(record['fasta_path'])
    _, fields_ = swell_from_bam(record['bam_path'], tiles, genomes, thresholds, min_pos, min_pos_total_zero)
    fields[0].extend(fields_[0])
    formatted_fields = [("%."+str(dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in fields[0]]
    # Add additional metadata given in the input table
    formatted_fields += [record.get(column) for column in metadata_headers]
    return "\t".join([str(x) for x in formatted_fields])


def swell_from_chunk(ichunk, genomes, metadata_headers, thresholds, dp, min_pos, min_pos_total_zero, clip):
    # Construct list of arguments to pass to executor
    args_list = [(record, genomes, metadata_headers, thresholds, dp, min_pos, min_pos_total_zero, clip) for record in ichunk]
    # Run swell as multiple processes
    # The context manager is exited once all processes are completed
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(swell_from_row, args_list)
    # 'results' contains an iterator of the process outputs, in the order that the processes were started
    return results


def swell_from_table(table_path, genomes, thresholds, dp, min_pos=None, min_pos_total_zero=False, clip=True):
    swell_headers = ["fasta_path", "header", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap", "bam_path", "num_pos", "mean_cov"]
    swell_headers += ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)]
    swell_headers += ["pc_tiles_medcov_gte%d" % x for x in sorted(thresholds)]
    swell_headers += ["tile_n", "tile_vector"]

    with open(table_path) as table:
        reader = csv.DictReader(table, delimiter='\t')
        table_headers = reader.fieldnames
        if table_headers:
            metadata_headers = [x for x in table_headers if not (x in set(swell_headers))]
            swell_headers.extend(metadata_headers)
            print("\t".join(swell_headers))

            # Iterate through the lines in the file, processing num_lines worth of lines at a time
            num_lines = 12
            ichunk = []
            for i, record in enumerate(reader):
                ichunk.append(record)
                if ((i + 1) % num_lines) == 0:
                    # Finished obtaining a chunk of lines from the file
                    # Now, process the chunk
                    ochunk = swell_from_chunk(ichunk, genomes, metadata_headers, thresholds, dp, min_pos, min_pos_total_zero, clip)
                    # The output is printed and we move to the next chunk
                    print('\n'.join(ochunk))
                    ichunk = []
            # Process the remaining partially filled chunk of lines
            if ichunk:
                ochunk = swell_from_chunk(ichunk, genomes, metadata_headers, thresholds, dp, min_pos, min_pos_total_zero, clip)
                print('\n'.join(ochunk))
                ichunk = []


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    fasta_parser = subparsers.add_parser("fasta")
    fasta_parser.add_argument("fasta_path")
    fasta_parser.add_argument("--summarise", action="store_true")
    fasta_parser.add_argument("--dp", default=2, type=int, required=False)
    fasta_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    bam_parser = subparsers.add_parser("bam")
    bam_parser.add_argument("bam_path")
    bam_parser.add_argument("--bed")
    bam_parser.add_argument("--ref", required=True, nargs='+')
    bam_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    bam_parser.add_argument("--min-pos", type=int, required=False)
    bam_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    bam_parser.add_argument("--no-clip", action="store_true")
    bam_parser.add_argument("--dp", default=2, type=int, required=False)
    bam_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    depth_parser = subparsers.add_parser("depth")
    depth_parser.add_argument("depth_path")
    depth_parser.add_argument("--bed")
    depth_parser.add_argument("--ref", required=True, nargs='+')
    depth_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    depth_parser.add_argument("--min-pos", type=int, required=False)
    depth_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    depth_parser.add_argument("--no-clip", action="store_true")
    depth_parser.add_argument("--dp", default=2, type=int, required=False)
    depth_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    fasta_bam_parser = subparsers.add_parser("fasta-bam")
    fasta_bam_parser.add_argument("fasta_path")
    fasta_bam_parser.add_argument("--summarise", action="store_true")
    fasta_bam_parser.add_argument("bam_path")
    fasta_bam_parser.add_argument("--bed")
    fasta_bam_parser.add_argument("--ref", required=True, nargs='+')
    fasta_bam_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    fasta_bam_parser.add_argument("--min-pos", type=int, required=False)
    fasta_bam_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    fasta_bam_parser.add_argument("--no-clip", action="store_true")
    fasta_bam_parser.add_argument("--dp", default=2, type=int, required=False)
    fasta_bam_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    fasta_depth_parser = subparsers.add_parser("fasta-depth")
    fasta_depth_parser.add_argument("fasta_path")
    fasta_depth_parser.add_argument("--summarise", action="store_true")
    fasta_depth_parser.add_argument("depth_path")
    fasta_depth_parser.add_argument("--bed")
    fasta_depth_parser.add_argument("--ref", required=True, nargs='+')
    fasta_depth_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    fasta_depth_parser.add_argument("--min-pos", type=int, required=False)
    fasta_depth_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    fasta_depth_parser.add_argument("--no-clip", action="store_true")
    fasta_depth_parser.add_argument("--dp", default=2, type=int, required=False)
    fasta_depth_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    table_parser = subparsers.add_parser("table")
    table_parser.add_argument("table_path")
    table_parser.add_argument("--ref", required=False, default=[], nargs='+')
    table_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    table_parser.add_argument("--min-pos", type=int, required=False)
    table_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    table_parser.add_argument("--no-clip", action="store_true")
    table_parser.add_argument("--dp", default=2, type=int, required=False)

    args = parser.parse_args()

    header = []
    fields = []

    if args.command:
        if (not args.command in ["fasta", "table"]) and args.bed:
            tiles = load_scheme(args.bed, not args.no_clip)
        else:
            tiles = {}

        if args.command == "fasta":
            if (not args.summarise):
                header_, fields_ = swell_from_fasta(args.fasta_path)
            else:
                header_, fields_ = summarise_swell_from_fasta(args.fasta_path)
            header.extend(header_)
            fields.extend(fields_)

        elif args.command == "bam":
            header_, fields_ = swell_from_bam(args.bam_path, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields.extend(fields_)
        
        elif args.command == "depth":
            header_, fields_ = swell_from_depth(args.depth_path, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields.extend(fields_)
        
        elif args.command == "fasta-bam":
            if (not args.summarise):
                header_, fields_ = swell_from_fasta(args.fasta_path)
            else:
                header_, fields_ = summarise_swell_from_fasta(args.fasta_path)  
            header.extend(header_)
            fields.extend(fields_)
            header_, fields_ = swell_from_bam(args.bam_path, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields[0].extend(fields_[0])
        
        elif args.command == "fasta-depth":
            if (not args.summarise):
                header_, fields_ = swell_from_fasta(args.fasta_path)
            else:
                header_, fields_ = summarise_swell_from_fasta(args.fasta_path)  
            header.extend(header_)
            fields.extend(fields_)
            header_, fields_ = swell_from_depth(args.depth_path, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields[0].extend(fields_[0])
        
        elif args.command == "table":
            swell_from_table(args.table_path, args.ref, args.thresholds, args.dp, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero, clip=not args.no_clip)
        
        if (not args.command == "table"):
            keys = []
            values = []
            if args.x:
                for key, value in args.x:
                    keys.append(key)
                    values.append(value)
            header.extend(keys)
            fields[0].extend(values)

            print("\t".join(header))
            for row in fields:
                row_s = [("%."+str(args.dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in row] # do not fucking @ me
                print("\t".join([str(x) for x in row_s]))


if __name__ == "__main__":
    main()