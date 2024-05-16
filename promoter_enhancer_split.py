# Code for assessing whether peaks are in the proximal or distal region
# Example usage:
# python promoter_enhancer_split.py arid1b_peaks_hg38.txt human ARID1B

import sys
import os
from matplotlib import pyplot as plt  # Basic plotting
import numpy as np
import seaborn as sns
import math

# Making a suitable gencode file
# awk -F '\t' '$3~/transcript/' gencode.vM23.annotation.gff3 > gencode.vM23.annotation_transcripts.txt
# awk -F '\t' '$3~/transcript/' gencode.v31.annotation.gtf > gencode.v31.annotation_transcripts.txt

def main():
    input_file = sys.argv[1]  # format chr \t start \t stop \t ...
    species = sys.argv[2]  # determines the promoter file to use
    tag_name = sys.argv[3]  # label used in PDFs
    gencode_promoter_file = ''
    
    # Note these files need to be ungzipped
    if species == 'human':
        gencode_promoter_file = 'gencode.v31.promoters.all.hg38.sort.txt'
    elif species == 'human_pc':
        gencode_promoter_file = 'gencode.v31.promoters.all.hg38.sort.pc.txt'
    elif species == 'mouse':
        gencode_promoter_file = 'gencode.vM23.promoters.all.mm10.txt'
    elif species == 'mouse_pc':
        gencode_promoter_file = 'gencode.vM23.promoters.all.mm10.pc.txt'

    fig_dir = ''

    # Runs the bedtools intersect command
    output_file = closest_tss(input_file, gencode_promoter_file)

    # Interprets the bedtools intersect file
    dist_to_tss(gencode_promoter_file, input_file, output_file, fig_dir, tag_name)


def count_file_lines(file_name_and_path):
    with open(file_name_and_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def closest_tss(input_file, gencode_file):
    # Use bedtools to find the closest TSS
    output_file = ''
    if input_file.endswith('.txt'):
        output_file = input_file.replace('.txt', '_gene.txt')
    else:
        output_file = f'{input_file}.gene'

    command = f'bedtools closest -a {input_file} -b {gencode_file} > {output_file}'

    print(command)
    os.system(command)

    return output_file


def dist_to_tss(gencode_promoter_file, input_file, output_file, fig_dir, tag_name):

    # Get the line length from the promoter file
    promoter_file = open(gencode_promoter_file, 'r')
    first_line = promoter_file.readline()
    tabs = first_line.rstrip().split('\t')
    prom_col_count = len(tabs)
    promoter_file.close()

    # Get the line length from the input file
    in_file = open(input_file, 'r')
    in_first_line = in_file.readline()
    in_tabs = in_first_line.rstrip().split('\t')
    in_col_count = len(in_tabs)
    in_file.close()

    # Calculate distance to TSS and assign promoter/enhancer
    include_chr = ['chrX', 'chrY', 'chrM']
    for i in range(22):
        include_chr.append(f'chr{i + 1}')

    # Work out file names and open files
    original_peak_count = count_file_lines(input_file)
    gene_file = open(output_file, 'r')

    dist_out_file = dist_nodup_out_file = ''
    if output_file.endswith('.txt'):
        dist_out_file = output_file.replace('.txt', '_dist.txt')
        dist_nodup_out_file = output_file.replace('.txt', '_dist_nodup.txt')
    else:
        dist_out_file = f'{output_file}.dist'
        dist_nodup_out_file = f'{output_file}.dist.nodup'

    dist = open(dist_out_file, 'w')
    nodup = open(dist_nodup_out_file, 'w')

    line_dict = {}
    dist_dict = {}
    type_dict = {}
    dir_dict = {}
    line_count = 0
    for line in gene_file:
        line_count += 1
        tab = line.rstrip().split('\t')

        # The first columns are from the input file
        chromo = tab[0]
        peak_start = int(tab[1])
        peak_end = int(tab[2])
        peak_key = f'{chromo}:{peak_start}-{peak_end}'

        # The next columns are from the promoter file
        prom_start = int(tab[in_col_count + 1])
        prom_end = int(tab[in_col_count + 2])
        strand = tab[in_col_count + 3]

        if chromo in include_chr:
            tss_pos = 0
            prom_pos = 0

            if strand == '+':
                tss_pos = prom_end
                prom_pos = prom_start
            elif strand == '-':
                tss_pos = prom_start
                prom_pos = prom_end
            else:
                print(f'ERR: Unrecognized strand: {strand}')

            peak_type = 'null'
            tss_dist = 99999999999
            direction = '.'

            # Peak overlaps the TSS
            if peak_start <= tss_pos and peak_end >= tss_pos:
                peak_type = 'promoter'
                tss_dist = 0

            # Peak is before the TSS
            elif peak_end < tss_pos:
                tss_dist = tss_pos - peak_end
                if peak_end < prom_pos:
                    peak_type = 'enhancer'
                    if strand == '+':
                        direction = 'upstream'
                    else:
                        direction = 'downstream'
                elif peak_end >= prom_pos:
                    peak_type = 'promoter'
                else:
                    print(f'ERR: Peak before TSS but relationship unclear')
                    print(f'Peak start: {peak_start}, Peak end: {peak_end}')
                    print(f'Prom start: {prom_start}, Prom end: {prom_end}')
                    print(f'Strand: {strand}, TSS pos: {tss_pos}')

            # Peak is after the TSS
            elif peak_start > tss_pos:
                tss_dist = peak_start - tss_pos
                if peak_start > prom_pos:
                    peak_type = 'enhancer'
                    if strand == '-':
                        direction = 'upstream'
                    else:
                        direction = 'downstream'
                elif peak_start <= prom_pos:
                    peak_type = 'promoter'
                else:
                    print(f'ERR: Peak after TSS but relationship unclear')
                    print(f'Peak start: {peak_start}, Peak end: {peak_end}')
                    print(f'Prom start: {prom_start}, Prom end: {prom_end}')
                    print(f'Strand: {strand}, TSS pos: {tss_pos}')

            # Other relationship
            else:
                print(f'ERR: Peak relationship unclear')
                print(f'Peak start: {peak_start}, Peak end: {peak_end}')
                print(f'Prom start: {prom_start}, Prom end: {prom_end}')
                print(f'Strand: {strand}, TSS pos: {tss_pos}')

            # Print out
            line_out = f'{line.rstrip()}\t{peak_type}\t{direction}\t{tss_dist}'
            if direction == 'downstream':
                line_out = f'{line.rstrip()}\t{peak_type}\t{direction}\t{-tss_dist}'

            dist.write(f'{line_out}\n')
            peak_dist_key = f'{peak_key}_{tss_dist}'
            line_dict[peak_dist_key] = line_out
            type_dict[peak_dist_key] = peak_type
            dir_dict[peak_dist_key] = direction
            if peak_key in dist_dict:
                dist_dict[peak_key].append(tss_dist)
            else:
                dist_dict[peak_key] = [tss_dist]

        else:
            tss_dist = '.'
            line_out = f'{line.rstrip()}\t.\t.\t.'
            dist.write(f'{line_out}\n')
            if peak_key in dist_dict:
                dist_dict[peak_key].append(tss_dist)
            else:
                dist_dict[peak_key] = [tss_dist]

            peak_dist_key = f'{peak_key}_{tss_dist}'
            line_dict[peak_dist_key] = line_out
            type_dict[peak_dist_key] = '.'
            dir_dict[peak_dist_key] = '.'

    # Print out without duplicates
    tss_dist_list = []
    prom_count = enhance_count = down_count = up_count = 0
    for peak_key in dist_dict:
        min_tss_dist = min(dist_dict[peak_key])
        peak_dist_key = f'{peak_key}_{min_tss_dist}'
        if dir_dict[peak_dist_key] == 'downstream':
            tss_dist_list.append(-min_tss_dist - 1)
        elif min_tss_dist != '.':
            tss_dist_list.append(min_tss_dist + 1)

        line_out = line_dict[peak_dist_key]
        nodup.write(f'{line_out}\n')
        if type_dict[peak_dist_key] == 'promoter':
            prom_count += 1
        elif type_dict[peak_dist_key] == 'enhancer':
            enhance_count += 1
            if dir_dict[peak_dist_key] == 'downstream':
                down_count += 1
            else:
                up_count += 1

    peak_count = len(dist_dict)
    missing = peak_count - prom_count - enhance_count

    prom_per = prom_count / peak_count
    enhance_per = enhance_count / peak_count
    down_per = down_count / peak_count
    up_per = up_count / peak_count
    prom_per_form = "{:.1%}".format(prom_per)
    enhance_per_form = "{:.1%}".format(enhance_per)
    down_per_form = "{:.1%}".format(down_per)
    up_per_form = "{:.1%}".format(up_per)

    print(f'{original_peak_count} peaks in the original file')
    print(f'Identified {peak_count} unique peaks with {prom_count} ({prom_per_form})',
          f'promoters and {enhance_count} ({enhance_per_form}) enhancers with {missing} difference in counts')

    # Plot the up, prom, down values
    bar_fig_out = ''
    out_file_simp = output_file.split('/')[-1]
    if out_file_simp.endswith('.txt'):
        hist_fig_out = f'{fig_dir}{out_file_simp.replace(".txt", "_hist.pdf")}'
        den_fig_out = f'{fig_dir}{out_file_simp.replace(".txt", "_density.pdf")}'
        lden_fig_out = f'{fig_dir}{out_file_simp.replace(".txt", "_log2density.pdf")}'
        bar_fig_out = f'{fig_dir}{out_file_simp.replace(".txt", "_barplot.pdf")}'
    else:
        hist_fig_out = f'{fig_dir}{out_file_simp}_hist.pdf'
        den_fig_out = f'{fig_dir}{out_file_simp}_density.pdf'
        lden_fig_out = f'{fig_dir}{out_file_simp}_log2density.pdf'
        bar_fig_out = f'{fig_dir}{out_file_simp}_barplot.pdf'

    # Make a barplot
    fig, ax = plt.subplots()
    plt.rcParams['pdf.use14corefonts'] = True
    ax.barh([tag_name], [up_per], 0.3, color='r')
    ax.barh([tag_name], [prom_per], 0.3, left=up_per, color='b')
    ax.barh([tag_name], [down_per], 0.3, left=up_per + prom_per, color='g')
    ax.set_xlabel(f'up: {up_per_form}, prom: {prom_per_form}, down: {down_per_form}')
    current_figure = plt.gcf()
    current_figure.savefig(bar_fig_out)
    plt.close()

    # Make a histogram
    fig, ax = plt.subplots()
    plt.rcParams['pdf.use14corefonts'] = True
    sns.set_style('whitegrid')
    sns.histplot(np.array(tss_dist_list), bins=20000)
    plt.xlim(5000, -5000)
    # xmin, xmax, ymin, ymax = plt.axis()
    # plt.vlines(0, ymin, ymax, colors='black')
    # plt.vlines(2000, ymin, ymax, colors='black')
    current_figure = plt.gcf()
    current_figure.savefig(hist_fig_out)
    plt.close()

    # Make a density plot
    fig, ax = plt.subplots()
    plt.rcParams['pdf.use14corefonts'] = True
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(tss_dist_list), bw_adjust=0.1)
    plt.xlim(100000, -100000)
    xmin, xmax, ymin, ymax = plt.axis()
    plt.vlines(0,ymin,ymax,colors='black')
    plt.vlines(2000,ymin,ymax,colors='black')
    current_figure = plt.gcf()
    current_figure.savefig(den_fig_out)
    plt.close()

    # Make a log-scaled density plot
    tss_dist_list_log = []
    for x in tss_dist_list:
        log_x = math.log(abs(x) + 1, 2)
        if x < 0:
            log_x = -log_x
        tss_dist_list_log.append(log_x)

    fig, ax = plt.subplots()
    plt.rcParams['pdf.use14corefonts'] = True
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(tss_dist_list_log))
    current_figure = plt.gcf()
    current_figure.savefig(lden_fig_out)
    plt.close()


if __name__ == '__main__':
    main()
