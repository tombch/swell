import sys
import numpy as np
import math
import argparse


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


def tile_depth_graph(tile_vector, show_values_on_right=False):
    bars_list, scale = scaled_bars(tile_vector)
    if show_values_on_right:
        depth_graph_string = f"TILE\tMEDIAN DEPTH\nNUM\t{scale}\n"
    else:
        depth_graph_string = f"TILE\tMEDIAN\nNUM\tDEPTH\t{scale}\n"
    if len(tile_vector) > 0:
        for i in range(len(tile_vector)):
            # tiles are sorted
            if show_values_on_right:
                depth_graph_string += f"{i+1}\t{bars_list[i]}{' ' if len(bars_list[i]) > 1 else ''}{int(round(tile_vector[i]))}\n"
            else:
                depth_graph_string += f"{i+1}\t{int(round(tile_vector[i]))}\t{bars_list[i]}\n"
    else:
        depth_graph_string += f"{None}\t{None}\n"
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--swell-output", nargs='+')
    parser.add_argument("--show-values-on-right", action='store_true', default=False)
    args = parser.parse_args()
    for output in args.swell_output:
        # Slightly hacky, could either do this or take the last argument of swell
        # Advantage for this way is that swell can still have graphs specified as arguments without interfering
        if ',' in output:
            tile_vector = np.fromstring(output, sep=",")
            break
    print(f"tile_vector = {tile_vector}\n")
    graph = tile_depth_graph(tile_vector, args.show_values_on_right)
    histogram = tile_histogram(tile_vector, args.show_values_on_right)
    print(f"{graph}\n\n{histogram}")


if __name__ == '__main__':
    main()