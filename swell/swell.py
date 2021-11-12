import argparse
import sys
import pysam
import numpy as np
import re
import math


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


def swell_from_fasta(fasta_path):
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
        from . import readfq # thanks heng
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

    return ["fasta_path", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid", "pc_ambiguous", "longest_gap", "longest_ungap"], [fasta_path, num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid, prop_ambiguous, max_gap, max_ungap]


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

    return tile_vector, ["bam_path", "num_pos", "mean_cov"] + ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)] + ["pc_tiles_medcov_gte%d" % x for x in sorted(thresholds)] + ["tile_n", "tile_vector"], [depth_path.replace(".depth", ""), n_positions, avg_cov] + threshold_counts_prop + tile_threshold_counts_prop + [len(tile_vector), tile_vector_str]


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


def scaled_bars(data):  
    max_val = int(max(data))
    marker = int('1' + '0'*(len(str(max_val))-1))
    space_increment = 10
    scale = "0" + " "*(space_increment - 1)
    for i in range(1, max_val):
        if i % marker == 0:
            i_str = str(i)
            i_len = len(i_str)
            scale += i_str + " "*(space_increment - i_len)
    bars_list = []
    for i in range(len(data)):
        bars_list.append(f"|{'='*int(space_increment*data[i]/marker)}")
    return bars_list, scale


def tile_depth_graph(tile_vector, tiles, show_values_on_right=False):
    bars_list, scale = scaled_bars(tile_vector)
    if show_values_on_right:
        depth_graph_string = f"TILE\tSTART\tEND\tMEDIAN DEPTH\nNUM\tPOS\tPOS\t{scale}\n"
    else:
        depth_graph_string = f"TILE\tSTART\tEND\tMEDIAN\nNUM\tPOS\tPOS\tDEPTH\t{scale}\n"
    # Tiles are already sorted
    for i in range(len(tile_vector)):
        if show_values_on_right:
            depth_graph_string += f"{i+1}\t{tiles[i][2]['inside_start']}\t{tiles[i][2]['inside_end']}\t{bars_list[i]}{' ' if len(bars_list[i]) > 1 else ''}{int(round(tile_vector[i]))}\n"
        else:
            depth_graph_string += f"{i+1}\t{tiles[i][2]['inside_start']}\t{tiles[i][2]['inside_end']}\t{int(round(tile_vector[i]))}\t{bars_list[i]}\n"
    return depth_graph_string[:-1]


def tile_histogram(tile_vector, show_values_on_right=False):
    depth_markers = [0, 10, 100, 1000, 10000, 100000, math.inf]
    hist_freqs = [0 for i in range(len(depth_markers) - 1)]
    for i in range(len(hist_freqs)):
        for depth in tile_vector:
            if depth_markers[i] <= depth < depth_markers[i + 1]:
                hist_freqs[i] += 1
    bars_list, scale = scaled_bars(hist_freqs)
    if show_values_on_right:
        histogram_string = f"FROM\tTO\tNUM TILES\nDEPTH\tDEPTH\t{scale}\n"
    else:
        histogram_string = f"FROM\tTO\tNUM\nDEPTH\tDEPTH\tTILES\t{scale}\n"
    for i in range(len(hist_freqs)):
        if show_values_on_right:
            histogram_string += f"{depth_markers[i]}\t{depth_markers[i+1]-1}\t{bars_list[i]}{' ' if len(bars_list[i]) > 1 else ''}{round(hist_freqs[i], 2)}\n"
        else:
            histogram_string += f"{depth_markers[i]}\t{depth_markers[i+1]-1}\t{round(hist_freqs[i], 2)}\t{bars_list[i]}\n"
    return histogram_string[:-1]


class ArgumentParserError(Exception):
    pass


class ErrorThrowingArgParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)


def main():
    parser = ErrorThrowingArgParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam", nargs=1, action='append') # bam file
    group.add_argument("--depth", nargs=1, action='append') # depth file produced from bam via samtools
    parser.add_argument("--bed", nargs=1, action='append', required=False)
    parser.add_argument("--fasta", nargs=1, action='append', required=False)
    parser.add_argument("--ref", required=True, nargs='+') # sequence name
    parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    parser.add_argument("--dp", default=2, type=int, required=False)
    parser.add_argument("--min-pos", type=int, required=False)
    parser.add_argument("--min-pos-allow-total-zero", action="store_true")
    parser.add_argument("-x", action="append", nargs=2, metavar=("key", "value",))
    parser.add_argument("--no-tile-clipping", action="store_true")
    parser.add_argument("--tile-depth-graph", action="store_true")
    parser.add_argument("--tile-histogram", action="store_true")
    args = parser.parse_args()

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

    fields = []
    header = []

    header_, fields_ = swell_from_fasta(args_fasta)
    header.extend(header_)
    fields.extend(fields_)

    if args_bam:
        tile_vector, header_, fields_ = swell_from_bam(args_bam, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
    elif args_depth:
        tile_vector, header_, fields_ = swell_from_depth(args_depth, tiles, args.ref, args.thresholds, min_pos=args.min_pos, min_pos_total_zero=args.min_pos_allow_total_zero)
    
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
    fields_s = [("%."+str(args.dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in fields] # do not fucking @ me
    print("\t".join([str(x) for x in fields_s]))

    graph = ""
    if args.tile_depth_graph:
        graph = tile_depth_graph(tile_vector, tiles)
        print(graph)
    
    histogram = ""
    if args.tile_histogram:
        histogram = tile_histogram(tile_vector)
        print(histogram)

if __name__ == "__main__":
    main()