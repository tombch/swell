import sys
import csv
import math
import random
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib.backends.backend_pdf import PdfPages


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
    # Tiles are already sorted
    for i in range(len(tile_vector)):
        if show_values_on_right:
            depth_graph_string += f"{i+1}\t{bars_list[i]}{' ' if len(bars_list[i]) > 1 else ''}{int(round(tile_vector[i]))}\n"
        else:
            depth_graph_string += f"{i+1}\t{int(round(tile_vector[i]))}\t{bars_list[i]}\n"
    depth_graph_title = "MEDIAN DEPTH PER TILE"
    depth_graph = depth_graph_string[:-1]
    return f"{'-' * len(depth_graph_title)}\n{depth_graph_title}\n{'-' * len(depth_graph_title)}\n{depth_graph}\n\n"


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
    histogram_title = "TILE HISTOGRAM"
    histogram = histogram_string[:-1]
    return f"{'-' * len(histogram_title)}\n{histogram_title}\n{'-' * len(histogram_title)}\n{histogram}"


def make_pdf(swell_data, pdf_path):
    fnt_size = 15
    mpl.style.use('ggplot')
    colour = "tab:orange"
    with PdfPages(pdf_path) as pdf:
        if "tile_vector" in swell_data[0].keys():
            tile_vector = swell_data[0]['tile_vector']
            bam_path = swell_data[0]['bam_path']
            f1 = plt.figure()
            ax_a = f1.add_subplot(111)
            ax_a.set(xlabel='tile', ylabel='median depth') 
            ax_a.bar([x+1 for x in range(len(tile_vector))], tile_vector, color=colour)
            f1.suptitle(f"bam_path = {bam_path}")
            f1.tight_layout()
            pdf.savefig() 
            plt.close()

            f2 = plt.figure()
            ax_b = f2.add_subplot(211)
            ax_b.set(xlabel='tile median depth', ylabel='# tiles')
            ax_b.hist(tile_vector, bins=[x for x in range(0, int(round(max(tile_vector) + 1)), 100)], color=colour)     
            ax_c = f2.add_subplot(212)
            ax_c.set(xlabel='tile median depth')
            ax_c.scatter(tile_vector, [random.random() for x in range(len(tile_vector))], s=10, color=colour)
            ax_c.get_yaxis().set_visible(False)
            f2.suptitle(f"bam_path = {bam_path}")
            f2.tight_layout()
            pdf.savefig() 
            plt.close()

        for row in swell_data:
            f3 = plt.figure()
            ax_d = f3.add_subplot(111)
            pie_data = {'pc_acgt' : None, 'pc_masked' : None, 'pc_invalid' : None, 'pc_ambiguous' : None}
            nonzero_pie_data = {}
            for k in pie_data.keys():
                pie_data[k] = row[k]
                if pie_data[k] > 0:
                    nonzero_pie_data[k] = pie_data[k]
            ax_d.pie(nonzero_pie_data.values(), colors=['tab:orange', 'lightgrey', 'black', 'tab:blue'], labels=nonzero_pie_data.keys(), startangle=90,  autopct='%1.1f%%')
            ax_d.title.set_text(f"fasta_path = {row['fasta_path']}\nbiosample_source_id = {row['biosample_source_id']}\nnum_seqs = {row['num_seqs']}")
            f3.tight_layout()
            pdf.savefig()
            plt.close()

        if len(swell_data) > 1:
            f4 = plt.figure()
            ax_e = f4.add_subplot(111)
            ax_e.set(xlabel='% acgt')
            birm_data = [swell_data[i]['pc_acgt'] for i in range(len(swell_data)) if swell_data[i]['sequencing_org'] == 'BIRM']
            ax_e.scatter(birm_data, [random.random() for x in range(len(birm_data))], s=10, color=colour)
            ax_e.get_yaxis().set_visible(False)
            f4.tight_layout()
            pdf.savefig() 
            plt.close()

            f5 = plt.figure()
            ax_f = f5.add_subplot(211)
            ax_f.set(xlabel='% acgt', ylabel='# samples')
            birm_data = [swell_data[i]['pc_acgt'] for i in range(len(swell_data)) if swell_data[i]['sequencing_org'] == 'BIRM']
            ax_f.hist(birm_data, bins=[x for x in range(0, 101, 5)], color=colour)     
            f5.tight_layout()
            pdf.savefig() 
            plt.close()

    print(f"PDF saved to: {pdf_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("swell_tsv")
    parser.add_argument("--show-values-on-right", action='store_true', default=False)
    parser.add_argument("--pdf-path")
    args = parser.parse_args()
    if args.swell_tsv == '-':
        f = sys.stdin.read().splitlines()
        swell_data = list(csv.DictReader(f, delimiter='\t'))
    else:
        with open(args.swell_tsv) as f:
            swell_data = list(csv.DictReader(f, delimiter='\t'))
    for row in swell_data:
        for k in row.keys():
            if k == 'tile_vector':
                row[k] = np.fromstring(row[k], sep=",")
            else:
                try:
                    row[k] = int(row[k])
                except ValueError:
                    try:
                        row[k] = float(row[k])
                    except ValueError:
                        pass
    if args.pdf_path:
        make_pdf(swell_data, args.pdf_path)


if __name__ == '__main__':
    main()