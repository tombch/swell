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

    scheme_fh.close()
    return new_tiles


def calculate_fasta_stats(seq):
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

    if seq:
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
    
    return [num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap]


def swell_from_fasta(fasta_path):
    '''
    Calculate fasta statistics given the path to a fasta/multifasta.
    '''
    if fasta_path == "-":
        fastas = readfq.readfq(sys.stdin)
    else:
        fastas = readfq.readfq(open(fasta_path))
        
    rows = []
    for name, seq, qual in fastas:
        rows.append([fasta_path, name] + calculate_fasta_stats(seq))

    if fasta_path != "-":
        fastas.close()
    return ["fasta_path", "header", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], rows


def swell_from_fasta_seq(seq, fasta_path="", header=""):
    '''
    Calculate fasta statistics directly from a sequence.
    '''
    rows = [[fasta_path, header] + calculate_fasta_stats(seq)]
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

    if fasta_path != "-":
        fastas.close()
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
    header, rows = swell_from_depth_iter(depth_fh, depth_path, tiles, genomes, thresholds, min_pos, min_pos_total_zero)
    depth_fh.close()
    return header, rows


def swell_from_bam(bam_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    # samtools 1.14
    # pysam 0.18.0 - much better than 0.16
    depth_iterable = (x.group(0)[:-1] for x in re.finditer('.*\n', pysam.depth('-a', bam_path))) # type: ignore
    return swell_from_depth_iter(depth_iterable, bam_path, tiles, genomes, thresholds, min_pos, min_pos_total_zero)


def swell_from_row(record, seq, header, genomes, metadata_header, thresholds, dp, min_pos, min_pos_total_zero, clip):
    # If the reference attached to the current record is not already in the provided list of references, add it
    ref = record.get('ref')
    if ref and (not ref in genomes):
        genomes.add(ref)

    # If a .bed file is provided, load the scheme for it (clipping by default)
    tiles = {}
    if record.get('bed_path'):
        tiles = load_scheme(record.get('bed_path'), clip)
    
    # Run swell fasta and swell bam on the paths given in the record. 
    if seq and header:
        # If the sequence and header were directly provided from a multifasta, use this instead of opening an individual fasta file
        _, fields = swell_from_fasta_seq(seq, record.get('fasta_path'), header)
    else:
        _, fields = swell_from_fasta(record.get('fasta_path'))
    _, fields_ = swell_from_bam(record.get('bam_path'), tiles, genomes, thresholds, min_pos, min_pos_total_zero)
    fields[0].extend(fields_[0])
    
    # Format the fields, and add additional metadata given in the input table
    formatted_fields = [("%."+str(dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in fields[0]]
    formatted_fields += [record.get(column) for column in metadata_header]
    # Then, return the outputted fields as a tab-separated string
    return "\t".join([str(x) for x in formatted_fields])


def fasta_fetch(fasta, header):
    if (not isinstance(fasta, type(None))) and (not isinstance(header, type(None))):
        return fasta.fetch(header)


def swell_from_table(table_path, max_workers, genomes, thresholds, dp, multi_fasta_path=None, min_pos=None, min_pos_total_zero=False, clip=True):
    swell_header = ["fasta_path", "header", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap", "bam_path", "num_pos", "mean_cov"]
    swell_header += ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)]
    swell_header += ["pc_tiles_medcov_gte%d" % x for x in sorted(thresholds)]
    swell_header += ["tile_n", "tile_vector"]

    # If a multifasta file is provided, fasta sequences will be pulled from here (instead of their individual fasta file)
    if multi_fasta_path:
        fasta = pysam.FastaFile(multi_fasta_path) # type: ignore
        # Dictionary storing fasta headers, with the pag names as keys
        pag_to_header = {ref.split('|')[0] : ref for ref in fasta.references} 
    else:
        fasta = None
        pag_to_header = {}

    with open(table_path) as table:
        reader = csv.DictReader(table, delimiter='\t')
        table_header = reader.fieldnames
        if table_header:
            # Attach all table headers that swell doesn't output, this is additional metadata
            metadata_header = [x for x in table_header if (not x in set(swell_header))]
            swell_header.extend(metadata_header)
            print("\t".join(swell_header))

            # Convert to a set for O(1) checks of genome existence
            genomes = set(genomes)
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                # 'submit' schedules the function to be executed
                # Returns a 'Future' object, which allows us to check if the process' result and/or if it is running/done
                results = [executor.submit(swell_from_row, record, fasta_fetch(fasta, pag_to_header.get(record.get('pag_name'))), pag_to_header.get(record.get('pag_name')), genomes, metadata_header, thresholds, dp, min_pos, min_pos_total_zero, clip) for record in reader]
                
                # Prints as processes are completed
                for f in concurrent.futures.as_completed(results):
                    print(f.result())


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
    table_parser.add_argument("--multi-fasta")
    table_parser.add_argument("--ref", required=False, default=[], nargs='+')
    table_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    table_parser.add_argument("--min-pos", type=int, required=False)
    table_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    table_parser.add_argument("--no-clip", action="store_true")
    table_parser.add_argument("--dp", default=2, type=int, required=False)
    table_parser.add_argument("--max-workers", default=24, type=int)

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
            swell_from_table(args.table_path, args.max_workers, args.ref, args.thresholds, args.dp, args.multi_fasta, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero, clip=not args.no_clip)
        
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