import argparse
import sys
import pysam
import numpy as np
import re
import math
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


def load_scheme(scheme_bed, no_clip=False):
    tiles_dict = {}
    scheme_fh = open(scheme_bed)
    for line in scheme_fh:
        ref, start, end, tile, pool = line.strip().split()
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
        ref, start, end, tile, pool = line.strip().split()
        scheme, tile, side = tile.split("_", 2)
        tile_tup = (scheme, tile, tiles_dict[tile])
        if tiles_dict[tile]["inside_start"] != -1 and tiles_dict[tile]["inside_end"] != -1 and tile not in tiles_seen:
            tiles_list.append(tile_tup)
            tiles_seen.add(tile)

    tiles_list = sorted(tiles_list, key=lambda x: int(x[1])) # sort by tile number
    # Clips by default
    if not no_clip:
        new_tiles = clip_tiles(tiles_list)
    else:
        new_tiles = tiles_list

    return new_tiles


def individual_swell_from_fasta(fasta_path):
    num_seqs = 0
    num_bases = 0
    num_acgt = 0
    num_masked = 0
    num_invalid = 0
    num_ambiguous = 0

    n_ungaps = []
    n_gaps = []
    curr_gap_len = 0
    curr_ungap_len = 0

    prop_acgt = 0
    prop_masked = 0
    prop_invalid = 0
    prop_ambiguous = 0
    max_gap = 0
    max_ungap = 0
    rows = []
    if fasta_path:
        if fasta_path == "-":
            heng_iter = readfq.readfq(sys.stdin)
        else:
            heng_iter = readfq.readfq(open(fasta_path))
        for name, seq, qual in heng_iter:
            num_seqs = 0
            num_bases = 0
            num_acgt = 0
            num_masked = 0
            num_invalid = 0
            num_ambiguous = 0

            n_ungaps = []
            n_gaps = []
            curr_gap_len = 0
            curr_ungap_len = 0

            prop_acgt = 0
            prop_masked = 0
            prop_invalid = 0
            prop_ambiguous = 0
            max_gap = 0
            max_ungap = 0

            num_seqs += 1
            for base in seq:
                num_bases += 1
                gap = 1

                if base.upper() in 'ACGT':
                    num_acgt += 1
                    gap = 0
                elif base.upper() in 'NX':
                    num_masked += 1
                elif base.upper() in 'WSMKRYBDHV':
                    num_ambiguous += 1
                    gap = 0
                else:
                    num_invalid += 1

                if gap:
                    if curr_ungap_len > 0:
                        n_ungaps.append(curr_ungap_len)
                        curr_ungap_len = 0
                    curr_gap_len += 1
                elif not gap:
                    if curr_gap_len > 0:
                        n_gaps.append(curr_gap_len)
                        curr_gap_len = 0
                    curr_ungap_len += 1

            if curr_gap_len > 0:
                n_gaps.append(curr_gap_len)
            elif curr_ungap_len > 0:
                n_gaps.append(curr_ungap_len)

            if num_bases > 0:
                prop_acgt = num_acgt / num_bases * 100.0
                prop_masked = num_masked / num_bases * 100.0
                prop_invalid = num_invalid / num_bases * 100.0
                prop_ambiguous = num_ambiguous / num_bases * 100.0

                if len(n_gaps) > 0:
                    max_gap = max(n_gaps)
                else:
                    max_gap = 0

                if len(n_ungaps) > 0:
                    max_ungap = max(n_ungaps)
                else:
                    max_ungap = 0
            else:
                prop_invalid = 100.0     
            # This will likely be replaced as I doubt it covers all cases
            name_split = name.split('/')
            biosample_source_id = name_split[1]
            sequencing_org = name_split[2].split(':')[0]
            rows.append([fasta_path, biosample_source_id, sequencing_org, num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap])

    return ["fasta_path", "biosample_source_id", "sequencing_org", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], rows


def group_swell_from_fasta(fasta_path):
    num_seqs = 0
    num_bases = 0
    num_acgt = 0
    num_masked = 0
    num_invalid = 0
    num_ambiguous = 0

    n_ungaps = []
    n_gaps = []
    curr_gap_len = 0
    curr_ungap_len = 0

    prop_acgt = 0
    prop_masked = 0
    prop_invalid = 0
    prop_ambiguous = 0
    max_gap = 0
    max_ungap = 0

    if fasta_path:
        if fasta_path == "-":
            heng_iter = readfq.readfq(sys.stdin)
        else:
            heng_iter = readfq.readfq(open(fasta_path))
        for name, seq, qual in heng_iter:
            num_seqs += 1
            for base in seq:
                num_bases += 1
                gap = 1

                if base.upper() in 'ACGT':
                    num_acgt += 1
                    gap = 0
                elif base.upper() in 'NX':
                    num_masked += 1
                elif base.upper() in 'WSMKRYBDHV':
                    num_ambiguous += 1
                    gap = 0
                else:
                    num_invalid += 1

                if gap:
                    if curr_ungap_len > 0:
                        n_ungaps.append(curr_ungap_len)
                        curr_ungap_len = 0
                    curr_gap_len += 1
                elif not gap:
                    if curr_gap_len > 0:
                        n_gaps.append(curr_gap_len)
                        curr_gap_len = 0
                    curr_ungap_len += 1

        if curr_gap_len > 0:
            n_gaps.append(curr_gap_len)
        elif curr_ungap_len > 0:
            n_gaps.append(curr_ungap_len)

        if num_bases > 0:
            prop_acgt = num_acgt / num_bases * 100.0
            prop_masked = num_masked / num_bases * 100.0
            prop_invalid = num_invalid / num_bases * 100.0
            prop_ambiguous = num_ambiguous / num_bases * 100.0

            if len(n_gaps) > 0:
                max_gap = max(n_gaps)
            else:
                max_gap = 0

            if len(n_ungaps) > 0:
                max_ungap = max(n_ungaps)
            else:
                max_ungap = 0
        else:
            prop_invalid = 100.0

    return ["fasta_path", "biosample_source_id", "sequencing_org", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], [fasta_path, "-", "-", num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap]


def swell_from_depth_iter(depth_iterable, depth_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    threshold_counters = {
        threshold: 0 for threshold in thresholds
    }
    tile_threshold_counters = {
        threshold: 0 for threshold in thresholds
    }
    n_positions = 0
    avg_cov = 0

    cursor = 0
    if tiles:
        tile_starts = [t[2]["inside_start"] for t in tiles] # dont use -1 for 1-pos depth files
        tile_ends = [t[2]["inside_end"] for t in tiles]

        closest_cursor = min(tile_starts)

        stat_tiles = [0 for t in tiles]
        tile_data = [[] for t in tiles]

    n_lines = 0
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

        if tiles:
            # Check for new open tiles
            if pos >= closest_cursor:
                for t_i, t_start in enumerate(tile_starts):
                    # If we are in the region of tile t_i and the tile is previously unvisited
                    if pos >= t_start and stat_tiles[t_i] == 0: 
                        # Change the state of tile t_i to open
                        stat_tiles[t_i] = 1 

                next_possible_min = []
                # Get the minimum tile start out of the remaining unvisited tiles 
                for t_i, t_start in enumerate(tile_starts):
                    if stat_tiles[t_i] == 0:
                        next_possible_min.append(t_start)
                try:
                    closest_cursor = min(next_possible_min)
                except ValueError:
                    closest_cursor = sys.maxsize

            # Handle open tiles
            for t_i, t_state in enumerate(stat_tiles):
                # If tile t_i is open, append the current coverage/depth value to t_i's list of coverages
                if t_state == 1:
                    tile_data[t_i].append(cov)
                # If our current position is at the end of t_i
                if tile_ends[t_i] <= pos:
                    # Mark tile t_i as closed
                    stat_tiles[t_i] = -1 

    tile_vector = []
    for t_i, (scheme_name, tile_num, tile) in enumerate(tiles):
        len_win = len(tile_data[t_i])
        mean_cov = np.mean(tile_data[t_i])
        median_cov = np.median(tile_data[t_i])
        
        tile_vector.append(median_cov)

        # Count tile means above threshold
        for threshold in threshold_counters:
            if median_cov >= threshold:
                tile_threshold_counters[threshold] += 1
        # print(depth_path, tile_num, tile[0], tile[1], scheme_name, mean_cov, median_cov, len_win)

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

    return ["bam_path", "num_pos", "mean_cov"] + ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)] + ["pc_tiles_medcov_gte%d" % x for x in sorted(thresholds)] + ["tile_n", "tile_vector"], [depth_path.replace(".depth", ""), n_positions, avg_cov] + threshold_counts_prop + tile_threshold_counts_prop + [len(tile_vector), tile_vector_str]


def swell_from_depth(depth_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    depth_fh = open(depth_path)
    return swell_from_depth_iter(depth_fh, depth_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False)


def swell_from_bam(bam_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False):
    depth_iterable = (x.group(0) for x in re.finditer('.*\n', pysam.depth('-a', bam_path)[:-1]))
    return swell_from_depth_iter(depth_iterable, bam_path, tiles, genomes, thresholds, min_pos=None, min_pos_total_zero=False)


# def swell_from_bam(bam_path, tiles, genome):
#     bam = pysam.AlignmentFile(bam_path)

#    for (scheme_name, tile_num, tile) in tiles:
#        tile_cover = bam.count_coverage(genome, tile[0]-1, tile[1],
#                quality_threshold=0, read_callback="all")
#        flat_tile_cover = np.array(tile_cover).sum(axis=0)

#        mean_cov = np.mean(flat_tile_cover)
#        median_cov = np.median(flat_tile_cover)
#        print(bam_path, tile_num, tile[0], tile[1], scheme_name, mean_cov, median_cov)


class ArgumentParserError(Exception):
    pass


class ErrorThrowingArgParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)


def main():
    parser = ErrorThrowingArgParser()
    subparsers = parser.add_subparsers(dest='command')

    group_parser = subparsers.add_parser("group")
    group_parser.add_argument("--fasta", nargs=1, action='append', required=False)
    group_exclusive = group_parser.add_mutually_exclusive_group(required=False)
    group_exclusive.add_argument("--bam", nargs=1, action='append') # bam file
    group_exclusive.add_argument("--depth", nargs=1, action='append') # depth file produced from bam via samtools
    group_parser.add_argument("--bed", nargs=1, action='append', required=False)
    group_parser.add_argument("--ref", required=False, nargs='+') # sequence name, also changed from required=True
    group_parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    group_parser.add_argument("--dp", default=2, type=int, required=False)
    group_parser.add_argument("--min-pos", type=int, required=False)
    group_parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    group_parser.add_argument("--no-tile-clipping", action="store_true")
    group_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    individual_parser = subparsers.add_parser("individual")
    individual_parser.add_argument("--fasta", nargs=1, action='append', required=False)
    individual_parser.add_argument("--dp", default=2, type=int, required=False)
    individual_parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))

    arguments = sys.argv[1:]
    sub_command = arguments[0]
    if sub_command != "group" and sub_command != "individual":
        # Default subcommand is to average samples in fasta
        arguments = ["group"] + arguments

    args = parser.parse_args(arguments)

    header = []
    fields = []

    if args.command == "group":
        called_once_args = {"--bam" : args.bam, "--depth" : args.depth, "--bed" : args.bed, "--fasta": args.fasta}
        for k in called_once_args.keys():
            if called_once_args[k] and len(called_once_args[k]) > 1:
                raise ArgumentParserError(f"{k} can only be specified once")
            elif called_once_args[k] and len(called_once_args[k]) == 1:
                called_once_args[k] = called_once_args[k][0][0]
        args_bam = called_once_args["--bam"]
        args_depth = called_once_args["--depth"]
        args_bed = called_once_args["--bed"]
        args_fasta = called_once_args["--fasta"]
        if args_bed:
            tiles = load_scheme(args_bed, args.no_tile_clipping)
        else:
            tiles = {}
        header_, fields_ = group_swell_from_fasta(args_fasta)
        header.extend(header_)
        fields.extend(fields_)
        if args_bam:
            header_, fields_ = swell_from_bam(args_bam, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields.extend(fields_)
        elif args_depth:
            header_, fields_ = swell_from_depth(args_depth, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
            header.extend(header_)
            fields.extend(fields_)
    elif args.command == "individual":
        called_once_args = {"--fasta": args.fasta}
        for k in called_once_args.keys():
            if called_once_args[k] and len(called_once_args[k]) > 1:
                raise ArgumentParserError(f"{k} can only be specified once")
            elif called_once_args[k] and len(called_once_args[k]) == 1:
                called_once_args[k] = called_once_args[k][0][0]
        args_fasta = called_once_args["--fasta"]
        header_, fields_ = individual_swell_from_fasta(args_fasta)
        header.extend(header_)
        fields.extend(fields_)

    keys = []
    values = []
    if args.x:
        for meta in args.x:
            keys.append(meta[0])
            values.append(meta[1])
    header.extend(keys)
    fields.extend(values)

    print("\t".join(header))
    if args.command == "group":
        fields_s = [("%."+str(args.dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in fields] # do not fucking @ me
        print("\t".join([str(x) for x in fields_s]))
    elif args.command == "individual":
        for row in fields:
            row_s = [("%."+str(args.dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in row] # do not fucking @ me
            print("\t".join([str(x) for x in row_s]))


if __name__ == "__main__":
    main()