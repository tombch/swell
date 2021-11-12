import sys
import numpy as np
import math


def scaled_bars(data):  
    max_val = int(max(data))
    marker = int('1' + '0'*(len(str(max_val))-1))
    space_increment = 10
    # scale = "0" + " "*(space_increment - 1)
    # for i in range(1, max_val):
    #     if i % marker == 0:
    #         i_str = str(i)
    #         i_len = len(i_str)
    #         scale += i_str + " "*(space_increment - i_len)
    bars_list = []
    for i in range(len(data)):
        bars_list.append(f"{'='*int(space_increment*data[i]/marker)}")
    return bars_list


def tile_depth_graph(tile_vector):
    bars_list = scaled_bars(tile_vector)
    depth_graph_string = f"TILE\tMEDIAN\nNUM\tDEPTH\n"
    if len(tile_vector) > 0:
        for i in range(len(tile_vector)):
            # tiles are sorted
            depth_graph_string += f"{i+1}\t{bars_list[i]}{' ' if len(bars_list[i]) > 0 else ''}{int(round(tile_vector[i]))}\n"
    else:
        depth_graph_string += f"{None}\t{None}\n"
    return depth_graph_string[:-1]


def tile_histogram(tile_vector):
    depth_markers = [0, 10, 100, 1000, 10000, 100000, math.inf]
    hist_freqs = [0 for i in range(len(depth_markers) - 1)]
    for i in range(len(hist_freqs)):
        for depth in tile_vector:
            if depth_markers[i] <= depth < depth_markers[i + 1]:
                hist_freqs[i] += 1
    bars_list = scaled_bars(hist_freqs)
    histogram_string = f"FROM\tTO\tNUM\nDEPTH\tDEPTH\tTILES\n"
    for i in range(len(hist_freqs)):
        histogram_string += f"{depth_markers[i]}\t{depth_markers[i+1]-1}\t{bars_list[i]}{' ' if len(bars_list[i]) > 0 else ''}{round(hist_freqs[i], 2)}\n"
    return histogram_string[:-1]


def main():
    for arg in sys.argv:
        # Slightly hacky, could either do this or take the last argument of swell
        # Advantage for this way is that swell can still have graphs specified as arguments without interfering
        if ',' in arg:
            tile_vector = np.fromstring(arg, sep=",")
            break
    print(f"tile vector: {tile_vector}\n")
    graph = tile_depth_graph(tile_vector)
    histogram = tile_histogram(tile_vector)
    print(f"{graph}\n\n{histogram}")


if __name__ == '__main__':
    main()