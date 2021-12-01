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


def get_scale_bound_and_marker(values):
    scale_bound = int('1' + '0'*(len(str(round(max(values))))))
    if max(values) < scale_bound / 2:
        scale_bound = int(scale_bound / 2)
    marker = int(scale_bound / 10)
    for i in range(scale_bound, 0, -marker):
        if max(values) + marker > i:
            break
        else:
            scale_bound -= marker            
    return scale_bound, marker


def plot_graph(y_values):
    graph_height, marker = get_scale_bound_and_marker(y_values)
    y_scale_val = (i for i in range(graph_height, -1, -marker))
    step_val = int(marker / 5)
    graph = []
    visited = [False for x in y_values]
    current_height_char = ":"
    for current_height in range(graph_height, 0, -step_val):
        current_line = []
        for i, x in enumerate(y_values):
            if not visited[i] and x >= current_height:
                current_line.append(current_height_char)
                visited[i] = True
            elif x >= current_height:
                current_line.append(':')
            # elif current_height % marker == 0:
            #     current_line.append('_')
            else:
                current_line.append(' ')
        if current_height % marker == 0:
            graph.append(current_line + ['|', '- ', str(next(y_scale_val))])
        else: 
            graph.append(current_line + ['|'])
        if current_height_char == ":":
            current_height_char = "."
        else:
            current_height_char = ":"
    graph.append(['=' for x in y_values] + ['|', '- ', str(next(y_scale_val))])
    return graph


def depth_x_axis(values):
    ticks = ""
    scale_string = f"{' ' * 4}"
    for x in values:
        if x % 5 == 0:
            scale_string += str(x)
            scale_string += f"{' ' * (5 - len(str(x)))}"
            ticks += "*"
        else:
            ticks += " "
    return f"{ticks}\n{scale_string}"


def depth_graph(x_values, y_values):
    graph = plot_graph(y_values)
    graph.append(depth_x_axis(x_values))
    graph_name = "MEDIAN DEPTH PER TILE".center(int(len(graph[-1])/2))
    bars = f"{len(graph_name.strip())*'-'}".center(int(len(graph[-1])/2))
    title = f"{bars}\n{graph_name}\n{bars}"
    graph.append(title)
    graph_string = '\n'.join([''.join(line) for line in graph])
    return graph_string


def get_counts(values, graph_length, step_val):
    bin_vals = [x for x in range(0, graph_length, step_val)]
    counts = [0 for i in range(len(bin_vals) - 1)]
    for i in range(len(counts)):
        for val in values:
            if bin_vals[i] <= val < bin_vals[i + 1]:
                counts[i] += 1
    return counts


def hist_x_axis(marker, graph_height, step_val):
    ticks = "<"
    scale_string = f""
    for x in range(0, graph_height, step_val):
        if x % marker == 0:
            scale_string += str(x)
            scale_string += f"{' ' * (10 - len(str(x)))}"
            ticks += f"{'-' * (10 - len('><'))}" + "><"
    return f"{ticks[:-1]}\n{scale_string}"


def histogram(x_values):
    graph_length, x_marker = get_scale_bound_and_marker(x_values)
    x_step_val = int(x_marker / 10)
    y_values = get_counts(x_values, graph_length, x_step_val)
    graph = plot_graph(y_values)
    graph.append(hist_x_axis(x_marker, graph_length, x_step_val))
    graph_name = "TILE HISTOGRAM".center(int(len(graph[-1])/2))
    bars = f"{len(graph_name.strip())*'-'}".center(int(len(graph[-1])/2))
    title = f"{bars}\n{graph_name}\n{bars}"
    graph.append(title)
    graph_string = '\n'.join([''.join(line) for line in graph])
    return graph_string


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
    parser.add_argument("--ascii-depth-graph", action="store_true")
    parser.add_argument("--ascii-histogram", action="store_true")
    parser.add_argument("--show-values-on-right", action="store_true", default=False)
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

    if args.ascii_depth_graph:
        tile_numbers = range(1, len(swell_data[0]['tile_vector'])+1)
        print(depth_graph(tile_numbers, swell_data[0]['tile_vector']))

    if args.ascii_histogram:
        print(histogram(swell_data[0]['tile_vector']))

    if args.pdf_path:
        make_pdf(swell_data, args.pdf_path)


if __name__ == '__main__':
    main()