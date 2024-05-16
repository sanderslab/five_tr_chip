# Code for assessing and visualizing the overlap between ChIP peak sets
# Example usage:
# python analysis_code_NeuroDevEpi.py -m 0.4 -t hg38_GW23 -g ARID1B,BCL11A,FOXP1,TBR1,TCF7L2 -p arid1b_peaks_hg38.txt,bcl11a_peaks_hg38.txt,foxp1_peaks_hg38.txt,tbr1_peaks_hg38.txt,tcf7l2_peaks_hg38.txt -s 3.9


# Code for assessing overlaps in Chip-seq Data
import matplotlib
from upsetplot import plot
from matplotlib import pyplot as plt  # Basic plotting
import pyupset as pyu
from upsetplot import from_memberships
from argparse import ArgumentParser
from colour import Color
import numpy as np
import matplotlib as mpl

# Directory with files
here = ''

# Directory for figures
fig_dir = ''

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-g', '--genes', type=str, required=True)  # comma seperated, e.g. ATAC,ARID1B,BCL11A
    parser.add_argument('-a', '--atac', type=str, required=False)  # ATAC-seq peak file
    parser.add_argument('-p', '--peakfiles', type=str, required=True)  # comma seperated, e.g. ATAC.bed,ARID1B.bed,BCL11A.bed
    parser.add_argument('-t', '--tag', type=str, required=True)  # e.g. human_hg38_GW23
    parser.add_argument('-m', '--max_per', type=float, required=False, default=0.65)  # e.g. 0.4 # maximum percentage for gradient plot
    parser.add_argument('-s', '--stat', type=float, required=False, default=3.9)  # e.g. 3.9 # log10(P) threshold for inclusion

    args = parser.parse_args()
    return args.genes, args.peakfiles, args.tag, args.atac, args.max_per, args.stat

# Main functions of the overlap code
def main(gene_list_str, peakfile_names, name_tag, atac_file, max_per, stat):

    file_count = 0

    peak_files = peakfile_names.split(',')
    gene_list = gene_list_str.split(',')
    for i in range(len(peak_files)):
        file_count += 1
        print(f'{file_count}) Gene: {gene_list[i]} with file: {peak_files[i]}')

    atac_mode = 0
    if atac_file:
        # name_tag = f'atac_{name_tag}'
        atac_mode = 1
        file_count += 1
        print(f'{file_count}) ATAC-seq peaks from file: {atac_file}')
        gene_list.append('ATAC')
        peak_files.append(atac_file)

    # Make the UCSC tracks
    make_ucsc(peak_files, gene_list, name_tag)

    # Human: Integrate all peaks into a castle plot for UCSC and analysis
    out_name = f'{here}castle_plot_peaks_{name_tag}_{stat}.txt'
    out = open(out_name, 'w')
    out.write('chr\tstart\tend\tunion_start\tunion_end\tunion_gene_count\tunion_genes\tpeak_count\tpeaks\tgene_count\tgenes\n')

    out_union_name = f'{here}castle_plot_union_peaks_{name_tag}_{stat}.txt'
    out_u = open(out_union_name, 'w')
    out_u.write('chr\tunion_start\tunion_end\tunion_gene_count\tunion_genes\n')

    bed_name = f'{here}castle_plot_peaks_{name_tag}_{stat}.bed'
    bed = open(bed_name, 'w')
    bed.write(f'track type=bedGraph name=\'Peak_Step\' description=\'{name_tag} Peak Step\'\n')

    key_metrics = f'{here}key_metrics_{name_tag}_{stat}.txt'
    keym = open(key_metrics, 'w')

    barplot = f'{fig_dir}bar_plot_peaks_{name_tag}_{stat}.pdf'

    # Add all the peaks to a single dictionary and sort list by chr/start
    pd, peaks_sort = merge_sort_peaks(peak_files, gene_list, barplot, stat)

    # Process the file to identify overlaps and print out
    gene_set_counts = process_merged_peaks(pd, peaks_sort, out, out_u, bed)

    out.close()
    bed.close()
    out_u.close()

    # Print out the number of peak unions per gene set
    data_list = []
    name_list = []
    total_union = 0
    multi_gene_union = 0
    four_plus_gene_union = 0
    for gene_set in sorted(gene_set_counts, key=gene_set_counts.get, reverse=True):
        # print(gene_set, gene_set_counts[gene_set])
        genes = gene_set.split(',')
        total_union += gene_set_counts[gene_set]
        if ',' in gene_set:
            multi_gene_union += gene_set_counts[gene_set]
            
        if gene_set.count(',') >= 3:
            four_plus_gene_union += gene_set_counts[gene_set]
            
        # Format the input for the upset plot, part 1
        name_list.append(genes)
        data_list.append(gene_set_counts[gene_set])

    # Make the table for the petal plot heatmap
    table_output_file_name = f'{here}table_for_heatmap_petal_{name_tag}_{stat}.txt'
    # factors = ['ARID1B', 'BCL11A', 'FOXP1', 'TBR1', 'TCF4']
    max_frac = table_for_heatmap_petal(gene_set_counts, gene_list, total_union, table_output_file_name)
        
    # Print out some key numbers to appear in the write up
    print(f'\nTotal union peaks: {total_union}')
    keym.write(f'Total union peaks\t{total_union}\n')
    percent_multi_gene_union = "{:.1%}".format(multi_gene_union / total_union)
    print(f'Total multigene union peaks: {multi_gene_union}, ({percent_multi_gene_union})')
    keym.write(f'Total multigene union peaks\t{multi_gene_union}\n')
    keym.write(f'Percent multigene union peaks\t{percent_multi_gene_union}\n')
    
    tbr_alone_count = 0
    if 'TBR1' in gene_set_counts:
        tbr_alone_count = gene_set_counts["TBR1"]

    percent_tbr1_union = "{:.1%}".format(tbr_alone_count / total_union)
    print(f'TBR1 only union peaks: {tbr_alone_count}, ({percent_tbr1_union})')
    keym.write(f'TBR1 only union peaks\t{tbr_alone_count}\n')
    keym.write(f'Percent TBR1 only union peaks\t{percent_tbr1_union}\n')

    non_tbr1_union = total_union - tbr_alone_count
    percent_multi_gene_non_tbr1_union = "{:.1%}".format(multi_gene_union / non_tbr1_union)
    print(f'Multigene non TBR1 union peaks: {multi_gene_union}, ({percent_multi_gene_non_tbr1_union})')
    keym.write(f'Multigene non TBR1 union peaks\t{multi_gene_union}\n')
    keym.write(f'Percent multigene non TBR1 union peaks\t{percent_multi_gene_non_tbr1_union}\n')

    percent_four_plus_gene_non_tbr1_union = "{:.1%}".format(four_plus_gene_union / non_tbr1_union)
    print(f'Four plus non TBR1 union peaks: {four_plus_gene_union}, ({percent_four_plus_gene_non_tbr1_union})')
    keym.write(f'Four plus non TBR1 union peaks\t{four_plus_gene_union}\n')
    keym.write(f'Percent four plus non TBR1 union peaks\t{percent_four_plus_gene_non_tbr1_union}\n')

    # Format the input for the upset plot, part 2
    upset_input = from_memberships(name_list, data=data_list)

    # plot(upset_input)
    plot(upset_input)
    current_figure = plt.gcf()
    fig_file_name = f'{fig_dir}upset_five_tfs_{name_tag}.pdf'
    current_figure.savefig(fig_file_name)

    # Petal plot of overlaps
    table_output_file_name = f'{here}table_for_heatmap_petal_{name_tag}.txt'
    petal_output_file_name = f'{fig_dir}heatmap_petal_{name_tag}_{stat}_max_exp.svg'
    petal(table_output_file_name, petal_output_file_name, max_per, "white", "red", gene_list)
    petal_output_file_name = f'{fig_dir}heatmap_petal_{name_tag}_{stat}_max_self.svg'
    petal(table_output_file_name, petal_output_file_name, max_frac, "white", "red", gene_list)

    # Assess observed number of five TF union peaks and four TF union peaks
    union = open(out_union_name, 'r')
    union_head = union.readline()
    union_4_5tf_name = f'{here}castle_plot_union_peaks_{name_tag}_{stat}_4_5tf.txt'
    union_4_5tf = open(union_4_5tf_name, 'w')

    union_4_5tf_bed_name = f'{here}castle_plot_union_peaks_{name_tag}_{stat}_4_5tf.bed'
    union_4_5tf_bed = open(union_4_5tf_bed_name, 'w')

    union_5tf_name = f'{here}castle_plot_union_peaks_{name_tag}_{stat}_5tf.txt'
    union_5tf = open(union_5tf_name, 'w')

    union_5tf_bed_name = f'{here}castle_plot_union_peaks_{name_tag}_{stat}_5tf.bed'
    union_5tf_bed = open(union_5tf_bed_name, 'w')

    observed_five_count = 0
    observed_four_count = 0
    observed_three_count = 0
    observed_two_count = 0
    observed_one_count = 0
    observed_none_count = 0
    total_peaks = 0
    for line in union:
        total_peaks += 1
        tab = line.rstrip().split()
        tf_count_in_peak = int(tab[3])

        if 'ATAC' in tab[4]:
            tf_count_in_peak = tf_count_in_peak - 1

        if tf_count_in_peak == 5:
            observed_five_count += 1
            union_4_5tf.write(line)
            union_5tf.write(line)
            first_three = "\t".join(tab[0:3])
            union_4_5tf_bed.write(f'{first_three}\n')
            union_5tf_bed.write(f'{first_three}\n')
        elif tf_count_in_peak == 4:
            observed_four_count += 1
            union_4_5tf.write(line)
            first_three = "\t".join(tab[0:3])
            union_4_5tf_bed.write(f'{first_three}\n')
        elif tf_count_in_peak == 3:
            observed_three_count += 1
        elif tf_count_in_peak == 2:
            observed_two_count += 1
        elif tf_count_in_peak == 1:
            observed_one_count += 1
        elif tf_count_in_peak == 0:
            observed_none_count += 1
            
    print(f'Found {total_peaks} total peaks')
    keym.write(f'Total peaks\t{total_peaks}\n')

    observed_five_count_per = "{:.1%}".format(observed_five_count / total_peaks)
    print(f'Found {observed_five_count} ({observed_five_count_per}) peaks with all five TFs binding')
    keym.write(f'Peaks with all five TFs\t{observed_five_count}\n')
    keym.write(f'Percent of peaks with all five TFs\t{observed_five_count_per}\n')

    observed_four_count_per = "{:.1%}".format(observed_four_count / total_peaks)
    print(f'Found {observed_four_count} ({observed_four_count_per}) peaks with four of five TFs binding')
    keym.write(f'Peaks with four of five TFs\t{observed_four_count}\n')
    keym.write(f'Percent of peaks with four of five TFs\t{observed_four_count_per}\n')

    observed_three_count_per = "{:.1%}".format(observed_three_count / total_peaks)
    print(f'Found {observed_three_count} ({observed_three_count_per}) peaks with three of five TFs binding')
    keym.write(f'Peaks with three of five TFs\t{observed_three_count}\n')
    keym.write(f'Percent of peaks with three of five TFs\t{observed_three_count_per}\n')

    observed_two_count_per = "{:.1%}".format(observed_two_count / total_peaks)
    print(f'Found {observed_two_count} ({observed_two_count_per}) peaks with two of five TFs binding')
    keym.write(f'Peaks with two of five TFs\t{observed_two_count}\n')
    keym.write(f'Percent of peaks with two of five TFs\t{observed_two_count_per}\n')

    observed_one_count_per = "{:.1%}".format(observed_one_count / total_peaks)
    print(f'Found {observed_one_count} ({observed_one_count_per}) peaks with one of five TFs binding')
    keym.write(f'Peaks with one of five TFs\t{observed_one_count}\n')
    keym.write(f'Percent of peaks with one of five TFs\t{observed_one_count_per}\n')

    observed_none_count_per = "{:.1%}".format(observed_none_count / total_peaks)
    print(f'Found {observed_none_count} ({observed_none_count_per}) peaks with none of five TFs binding')
    keym.write(f'Peaks with none of five TFs\t{observed_none_count}\n')
    keym.write(f'Percent of peaks with none of five TFs\t{observed_none_count_per}\n')

    observed_four_five_count_per = "{:.1%}".format((observed_four_count + observed_five_count) / total_peaks)
    print(f'Found {observed_four_count + observed_five_count} ({observed_four_five_count_per}) peaks with four of more TFs binding')
    keym.write(f'Peaks with at least four TFs\t{observed_five_count}\n')
    keym.write(f'Percent of peaks with at least four TFs\t{observed_four_five_count_per}\n')

    union_4_5tf.close()
    union_4_5tf_bed.close()
    keym.close()


# key subroutines for assessing peak overlaps

# Make a UCSC track 
def make_ucsc(peak_files, gene_list, name_tag):
    if len(peak_files) != len(gene_list):
        print(f'Error: there are {len(peak_files)} and {len(gene_list)} genes, these numbers should match')
    
    for file_num in range(len(peak_files)):
        
        peak_file = peak_files[file_num]
        peaks = open(peak_file, 'r')

        the_gene = gene_list[file_num]
        
        peak_file_bed = f'{peak_file}_{name_tag}.bed'
        out = open(peak_file_bed, 'w')
        out.write(f'track name=\'{the_gene}_Peaks\' description=\'{name_tag} {the_gene} Peaks\'\n')

        line_count = 0
        for line in peaks:
            line_count += 1
            tab = line.rstrip().split('\t')
            peak_name = f'{name_tag}_{line_count}'
            if len(tab) >= 4:
                peak_name = tab[3].split('_')[-1]
            first_three_col = '\t'.join(tab[:3])
            outline = f'{first_three_col}\t{peak_name}\n'
            out.write(outline)

    out.close()
    peaks.close()

# Add the position to the coord dict
def add_coord(coord, peak, pos, place):
    if pos in coord:
        coord[pos] = f'{coord[pos]}|{place}_{peak}'
    else:
        coord[pos] = f'{place}_{peak}'
        
    return coord
    
# Process overlapping peaks
def process_last_matches(coord, chromo, gene_set, pdict):
    out_list = []
    union_list = []
    bed_list = []
    gene_here = {}
    peak_here = {}
    max_gene_count = len(gene_set)
    max_gene_list_str = ','.join(sorted(list(gene_set)))
    
    if max_gene_count >= 1:
        coord_sort = sorted(list(coord))
        union_start = coord_sort[0]
        union_end = coord_sort[-1]
        for i in range(len(coord_sort)-1):
            pos = coord_sort[i]
            next_pos = coord_sort[i+1]
            features = coord[pos].split('|')
            for feature in features:
                if feature.startswith('s_'):
                    peak = feature.replace('s_', '')
                    peak_here[peak] = 5
                    gene_here[pdict[peak]['gene']] = 5
                elif feature.startswith('e_'):
                    peak = feature.replace('e_', '')
                    del peak_here[peak]
                    gene_here.pop(pdict[peak]['gene'], None)

            peak_list_str = ','.join(sorted(list(peak_here)))
            gene_list_str = ','.join(sorted(list(gene_here)))
            gene_count = len(gene_here)
            peak_count = len(peak_here)
            outline = f'{chromo}\t{pos}\t{next_pos}\t{union_start}\t{union_end}\t{max_gene_count}\t{max_gene_list_str}\t{peak_count}\t{peak_list_str}\t{gene_count}\t{gene_list_str}\n'
            bedline = f'{chromo}\t{pos}\t{next_pos}\t{gene_count}\n'
            out_list.append(outline)
            bed_list.append(bedline)

        outline_union = f'{chromo}\t{union_start}\t{union_end}\t{max_gene_count}\t{max_gene_list_str}\n'
        union_list.append(outline_union)
        
    return out_list, bed_list, union_list

# Add all the peaks to a single dictionary and sort list by chr/start
def merge_sort_peaks(peak_files, gene_list, fig_file_name, stat):
    pd = {}
    uni_peak = {}
    peak_names = []
    peak_counts = []
    for i in range(len(peak_files)):
        peak_file = peak_files[i]
        the_gene = gene_list[i]
        peak_count = 0

        # Get the gene name
        # the_gene = peak_file.split('/')[-1].split('_')[name_split_factor].upper()

        peaks = open(peak_file, 'r')
        for line in peaks:
            
            tab = line.rstrip().split('\t')
            peak_name = f'{the_gene}_{peak_count}'
#             print(peak_name)
            if the_gene == 'ATAC':
                peak_count += 1
                pd[peak_name] = {'gene': the_gene, 'chr': tab[0], 'start': int(tab[1]), 'end': int(tab[2])}
            elif len(tab) >= 8:
                peak_stat = float(tab[7])

                # Add to dictionary
                if peak_stat >= stat:
                    peak_count += 1
                    pd[peak_name] = {'gene': the_gene, 'chr': tab[0], 'start': int(tab[1]), 'end': int(tab[2])}

        print(f'Loaded {peak_count} peaks for the gene {the_gene} at logP over {stat}')
        peak_names.append(the_gene)
        peak_counts.append(peak_count)

    # Sort peaks by start position then chromosome
    peaks = list(pd.keys())
    peaks_sort = sorted(peaks, key=lambda x: (pd[x]['chr'], pd[x]['start']))

    # Make a barplot
    fig, ax = plt.subplots()
    plt.rcParams['pdf.use14corefonts'] = True
    col=('#377EB8', '#FF7F00', '#4DAF4A', '#984EA3', '#E41A1C')
    peak_num = ax.bar(peak_names, peak_counts, 0.3, color=col)
    autolabel(ax, peak_num)
    ax.set_ylabel('Number of peaks')
    ax.set_ylim([0, 105000])
    current_figure = plt.gcf()
    current_figure.savefig(fig_file_name)
    
    return pd, peaks_sort


def autolabel(ax, rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1000, f'{height}',
                ha='center', va='bottom', rotation='vertical')


# Process the file to identify overlaps and print out
def process_merged_peaks(pdict, peaks_sort, out, out_u, bed):
    coord = {}
    gene_set = {}
    gene_set_counts = {}
    max_end = 0
    last_chromo = '.'
    for peak in peaks_sort:
        chromo = pdict[peak]['chr']
        start = int(pdict[peak]['start'])
        end = int(pdict[peak]['end'])
        gene = pdict[peak]['gene']

        if chromo == last_chromo and start <= max_end:
            # Match!
            if end > max_end:
                max_end = end

            # Add to coord
            coord = add_coord(coord, peak, start, 's')
            coord = add_coord(coord, peak, end, 'e')
            if gene not in gene_set:
                gene_set[gene] = 5

        else:
            # No match
            out_list, bed_list, union_list = process_last_matches(coord, last_chromo, gene_set, pdict)
            genes = ','.join(sorted(list(gene_set)))
            if genes in gene_set_counts:
                gene_set_counts[genes] += 1
            elif genes != '':
                gene_set_counts[genes] = 1

            for line_o in out_list:
                out.write(line_o)

            for line_b in bed_list:
                bed.write(line_b)

            for line_u in union_list:
                out_u.write(line_u)

            # Add to coord
            coord = {}  # clear coord
            coord = add_coord(coord, peak, start, 's')
            coord = add_coord(coord, peak, end, 'e')
            last_chromo = chromo
            max_end = end
            gene_set = {gene: 5}

    # Process last match
    out_list, bed_list, union_list = process_last_matches(coord, last_chromo, gene_set, pdict)
    genes = ','.join(sorted(list(gene_set)))
    if genes in gene_set_counts:
        gene_set_counts[genes] += 1
    elif genes != '':
        gene_set_counts[genes] = 1

    for line_o in out_list:
        out.write(line_o)

    for line_b in bed_list:
        bed.write(line_b)

    for line_u in union_list:
        out_u.write(line_u)
        
    return gene_set_counts


def count_file_lines(file_name_and_path):
    with open(file_name_and_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Make a table for a petal plot heatmap
def table_for_heatmap_petal(gene_set_counts, factors, total, table_output_file_name):
    table_out = open(table_output_file_name, 'w')
    
    max_frac = 0
    for gene_set in gene_set_counts:
        
        # Get the code
        code = ['0'] * len(factors)
        for factor_count in range(len(factors)):

            if factors[factor_count] in gene_set:
                code[factor_count] = '1'
                
        code_str = ''.join(code)
        
        # Get the count and fraction
        count = gene_set_counts[gene_set]
        fraction = count / total
        if fraction > max_frac:
            max_frac = fraction
        
        table_out.write(f'{gene_set}\t{code_str}\t{count}\t{fraction}\n')

    
    table_out.close()
    return max_frac

# Select colors for images
def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)

    # Prevent values out of range from exceeding the color scale
    if mix > 1:
        mix = 1

    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

# Make the petal plot
def petal(infile_name, outfile_name, max_per, c1, c2, gene_list):
    
    color_resolution = 1000
    max_factor = 1 / max_per  # e.g. 0.406

    out = open(outfile_name, 'w')
    table = open(infile_name, 'r')

    out.write('<?xml version="1.0" encoding="utf-8"?>\n')
    out.write('<!-- Generator: Adobe Illustrator 25.1.0, SVG Export Plug-In . SVG Version: 6.00 Build 0)  -->\n')
    out.write('<svg version="1.1" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px"\n')
    out.write('  viewBox="0 0 576 576" style="enable-background:new 0 0 576 576;" xml:space="preserve">\n')
    out.write('<style type="text/css">\n')
    out.write('\t.st1{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')

    # Initialize code dictionary
    combo = ['10000', '01000', '00100', '00010', '00001', '11000', '10100', '10010', '10001', '01100', '01010', '01001', '00110', 
    '00101', '00011', '11100', '11010', '11001', '10110', '10101', '10011', '01110', '01101', '01011', '00111', '11110', '11101',
    '11011', '10111', '01111', '11111']
    code_to_count = {}
    code_to_color = {}
    for code in combo:
        code_to_count[code] = 0
        code_to_color[code] = c1

    # Modify values from table
    count_text_out = []
    for line in table:
        tab = line.rstrip().split('\t')
        code = tab[1]
        code_to_count[code] = int(tab[2])
        
        color_out = c1
        color_frac = float(tab[3]) * max_factor
        color_out = colorFader(c1,c2,color_frac)
        # print(color_frac)
        code_to_color[code] = color_out

    # Color and counts
    for code in combo:
        color_code = code_to_color[code]
        count = code_to_count[code]
        if code == '10000':
            out.write('\t.st2{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 120.5335 248.7162)" class="st34 st35">{count}</text>')
        elif code == '01000':
            out.write('\t.st30{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 278.6256 133.2766)" class="st34 st35">{count}</text>')
        elif code == '00100':
            out.write('\t.st31{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 441.5117 248.7159)" class="st34 st35">{count}</text>')
        elif code == '00010':
            out.write('\t.st18{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 369.7153 447.2975)" class="st34 st35">{count}</text>')
        elif code == '00001':
            out.write('\t.st15{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 187.6006 447.2975)" class="st34 st35">{count}</text>')
        elif code == '11000':
            out.write('\t.st6{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 186.6005 212.9799)" class="st34 st35">{count}</text>')
        elif code == '10100':
            out.write('\t.st3{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 148.0147 303.1999)" class="st34 st35">{count}</text>')
        elif code == '10010':
            out.write('\t.st19{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 381.3955 389.2256)" class="st34 st35">{count}</text>')
        elif code == '10001':
            out.write('\t.st10{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 161.6805 363.3396)" class="st34 st35">{count}</text>')
        elif code == '01100':
            out.write('\t.st28{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 344.9078 185.8387)" class="st34 st35">{count}</text>')
        elif code == '01010':
            out.write('\t.st29{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 246.6725 179.452)" class="st34 st35">{count}</text>')
        elif code == '01001':
            out.write('\t.st12{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 232.0956 424.2198)" class="st34 st35">{count}</text>')
        elif code == '00110':
            out.write('\t.st22{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 412.5907 323.2693)" class="st34 st35">{count}</text>')
        elif code == '00101':
            out.write('\t.st26{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 385.5107 225.0196)" class="st34 st35">{count}</text>')
        elif code == '00011':
            out.write('\t.st16{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 294.0507 438.8972)" class="st34 st35">{count}</text>')
        elif code == '11100':
            out.write('\t.st4{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 179.8049 259.4544)" class="st34 st35">{count}</text>')
        elif code == '11010':
            out.write('\t.st5{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 211.1205 202.4147)" class="st34 st35">{count}</text>')
        elif code == '11001':
            out.write('\t.st11{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 211.1205 386.3497)" class="st34 st35">{count}</text>')
        elif code == '10110':
            out.write('\t.st21{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 399.7697 342.9097)" class="st34 st35">{count}</text>')
        elif code == '10101':
            out.write('\t.st8{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 160.6805 342.9097)" class="st34 st35">{count}</text>')
        elif code == '10011':
            out.write('\t.st17{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 339.5356 396.0097)" class="st34 st35">{count}</text>')
        elif code == '01110':
            out.write('\t.st25{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 284.6205 195.9897)" class="st34 st35">{count}</text>')
        elif code == '01101':
            out.write('\t.st27{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 355.9945 204.5796)" class="st34 st35">{count}</text>')
        elif code == '01011':
            out.write('\t.st14{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 280.6556 427.3506)" class="st34 st35">{count}</text>')
        elif code == '00111':
            out.write('\t.st23{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 381.3955 267.8547)" class="st34 st35">{count}</text>')
        elif code == '11110':
            out.write('\t.st0{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 221.8954 226.0797)" class="st34 st35">{count}</text>')
        elif code == '11101':
            out.write('\t.st9{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 187.6006 331.6697)" class="st34 st35">{count}</text>')
        elif code == '11011':
            out.write('\t.st13{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 284.6205 402.5913)" class="st34 st35">{count}</text>')
        elif code == '10111':
            out.write('\t.st20{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 377.3163 325.96)" class="st34 st35">{count}</text>')
        elif code == '01111':
            out.write('\t.st24{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 331.438 220.9398)" class="st34 st35">{count}</text>')
        elif code == '11111':
            out.write('\t.st7{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 278.998 303.2002)" class="st34 st35">{count}</text>')
        

    out.write('\t.st32{font-family:\'Helvetica-Oblique\';}\n')
    out.write('\t.st33{font-size:24px;}\n')
    out.write('\t.st34{font-family:\'Helvetica\';}\n')
    out.write('\t.st35{font-size:12px;}\n')
    out.write('\t.st36{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
    out.write('\t.st37{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
    out.write('</style>\n')

    # Locations
    out.write('<path class="st0" d="M209.86,275.16c-1.98-12.51-2.81-24.69-2.43-36.25c0.31-9.35,1.41-18.11,3.24-26.18l0,0\n')
    out.write(' c10.53-6.27,22.03-11.89,34.24-16.74l0,0c15.59,3.14,31.49,7.86,47.19,14.05l0,0c-25.18,12.73-50.16,31.37-72.16,54.15\n')
    out.write(' C216.48,267.78,213.12,271.44,209.86,275.16L209.86,275.16z"/>\n')
    out.write('<rect class="st1" width="576" height="576"/>\n')
    out.write('<path class="st2" d="M161.04,355.52c-13.44-11.52-24.99-23.96-34.06-36.85c-21.4-30.43-27.19-60.87-16.12-84.63\n')
    out.write(' c11.08-23.76,38.13-38.88,75.19-42.05c5.12-0.44,10.38-0.64,15.74-0.62l0,0c-4.33,11.02-8.03,22.66-11,34.71l0,0\n')
    out.write(' c-1.56,1.21-3.08,2.44-4.56,3.68c-28.54,23.85-42.02,51.76-37.47,77.57c2.18,12.35,8.39,23.73,18.12,33.64l0,0\n')
    out.write(' C164.72,345.83,162.77,350.69,161.04,355.52L161.04,355.52z"/>\n')
    out.write('<path class="st3" d="M190.94,226.02c-2.39,9.68-4.31,19.62-5.73,29.73c-2.66,18.9-3.46,37.72-2.48,55.76l0,0\n')
    out.write(' c-6.05,9.71-11.32,19.56-15.7,29.4l0,0c-9.72-9.9-15.94-21.29-18.12-33.64c-4.55-25.81,8.93-53.72,37.47-77.57\n')
    out.write(' C187.86,228.46,189.38,227.23,190.94,226.02L190.94,226.02z"/>\n')
    out.write('<path class="st4" d="M210.64,212.98c-1.84,8.07-2.94,16.83-3.24,26.19c-0.38,11.57,0.45,23.75,2.43,36.25l0,0\n')
    out.write(' c-10.2,11.67-19.29,23.93-27.09,36.43l0,0c-0.98-18.04-0.17-36.86,2.48-55.76c1.42-10.11,3.34-20.05,5.73-29.73l0,0\n')
    out.write(' C197.03,221.63,203.62,217.16,210.64,212.98L210.64,212.98z"/>\n')
    out.write('<path class="st5" d="M210.61,212.73c1.71-7.5,4.06-14.41,7.01-20.63l0,0c8.94,0.73,18.09,2.05,27.36,3.92l0.02-0.08\n')
    out.write(' C232.73,200.79,221.18,206.43,210.61,212.73L210.61,212.73z"/>\n')
    out.write('<path class="st6" d="M210.62,212.7c-7.02,4.18-13.6,8.65-19.68,13.38l0,0c2.97-12.05,6.67-23.69,11-34.71l0,0\n')
    out.write(' c5.14,0.02,10.38,0.25,15.69,0.69l0,0C214.68,198.29,212.33,205.19,210.62,212.7L210.62,212.7z"/>\n')
    out.write('<path class="st7" d="M382.83,267.21c-0.46,5.35-1.06,10.73-1.82,16.11c-4.49,31.97-14.02,62.29-27.31,87.69l0,0\n')
    out.write(' c-5.99,1.45-12.07,2.72-18.22,3.8c-29.93,5.28-59.85,5.86-86.84,1.91l0,0c-3.39-5.47-6.62-11.11-9.69-16.87\n')
    out.write(' c-14.75-27.74-24.72-57.02-29.1-84.71l0,0c3.25-3.72,6.62-7.39,10.08-10.98c22-22.78,46.98-41.42,72.16-54.15l0,0\n')
    out.write(' c5.56,2.2,11.1,4.58,16.6,7.14C337.13,230.42,362.69,247.77,382.83,267.21L382.83,267.21z"/>\n')
    out.write('<path class="st8" d="M167.04,340.98c4.38-9.84,9.65-19.69,15.7-29.4l0,0c0.88,16.27,3.22,31.91,6.95,46.36l0,0\n')
    out.write(' C180.82,352.97,173.22,347.28,167.04,340.98L167.04,340.98z"/>\n')
    out.write('<path class="st9" d="M209.99,275.17c4.38,27.69,14.35,56.97,29.1,84.71c3.07,5.76,6.3,11.4,9.69,16.87l0,0\n')
    out.write(' c-13.51-1.97-26.28-5.08-37.95-9.3c-7.61-2.75-14.62-5.93-20.97-9.48l0,0c-3.73-14.46-6.07-30.09-6.95-46.36l0,0\n')
    out.write(' C190.7,299.1,199.79,286.84,209.99,275.17L209.99,275.17z"/>\n')
    out.write('<path class="st10" d="M166.93,340.91c-2.17,4.86-4.12,9.72-5.84,14.56c10.96,9.39,23.17,18.17,36.32,26.09\n')
    out.write(' c-0.19-0.45-0.37-0.91-0.56-1.36c-2.83-7.04-5.26-14.51-7.28-22.33l0,0C180.71,352.9,173.11,347.21,166.93,340.91L166.93,340.91z"/>\n')
    out.write('<path class="st11" d="M248.61,376.73c8.52,13.77,18,26.56,28.13,38.02c-17.22-4.33-34.65-10.52-51.62-18.43\n')
    out.write(' c-9.56-4.46-18.79-9.37-27.6-14.68c-0.19-0.45-0.37-0.91-0.56-1.36c-2.83-7.04-5.26-14.51-7.28-22.33l0,0\n')
    out.write(' c6.35,3.55,13.36,6.73,20.97,9.48C222.33,371.65,235.11,374.75,248.61,376.73L248.61,376.73z"/>\n')
    out.write('<path class="st12" d="M259.06,439.81c-25.61-3.6-47.74-24.53-61.64-58.24c8.81,5.31,18.04,10.22,27.6,14.68\n')
    out.write(' c16.97,7.91,34.4,14.11,51.62,18.43c6.02,6.81,12.27,13.15,18.69,18.94C283.32,439.32,271.03,441.49,259.06,439.81z"/>\n')
    out.write('<path class="st13" d="M314.85,421.14c-12.34-1.1-25.08-3.29-37.94-6.53c-10.13-11.46-19.61-24.24-28.13-38.02l0,0\n')
    out.write(' c26.98,3.94,56.9,3.36,86.83-1.91c6.15-1.08,12.23-2.35,18.22-3.8c-5.48,10.48-11.61,20.13-18.28,28.71\n')
    out.write(' C328.98,408.03,322.03,415.24,314.85,421.14L314.85,421.14z"/>\n')
    out.write('<path class="st14" d="M314.68,421.28c-6.27,5.15-12.72,9.31-19.26,12.41c-6.42-5.79-12.67-12.13-18.69-18.94\n')
    out.write(' C289.6,417.98,302.34,420.17,314.68,421.28L314.68,421.28z"/>\n')
    out.write('<path class="st15" d="M308.75,444.71c-17.25,10.81-34.94,18.99-52.13,23.96c-35.74,10.32-66.41,5.89-85.26-12.31\n')
    out.write(' c-18.86-18.21-24.35-48.71-15.28-84.78c1.35-5.36,3-10.77,4.94-16.22c10.96,9.39,23.17,18.17,36.32,26.09\n')
    out.write(' c13.9,33.72,36.03,54.64,61.64,58.24c11.97,1.68,24.26-0.48,36.26-6.19C299.67,437.51,304.18,441.25,308.75,444.71L308.75,444.71z"\n')
    out.write(' />\n')
    out.write('<path class="st16" d="M314.57,421.14c8.76,0.78,17.33,1.01,25.59,0.66c-10.13,8.57-20.65,16.27-31.34,22.96l0,0\n')
    out.write(' c-4.57-3.46-9.08-7.2-13.51-11.2C301.85,430.45,308.3,426.3,314.57,421.14L314.57,421.14z"/>\n')
    out.write('<path class="st17" d="M314.69,421.21c7.18-5.9,14.13-13.11,20.7-21.55c6.67-8.58,12.8-18.22,18.28-28.71\n')
    out.write(' c15.44-3.73,30.26-8.63,44.08-14.52c-9.98,15.64-22,30.86-35.65,44.99c-7,7.25-14.31,14.09-21.83,20.44\n')
    out.write(' C332.01,422.21,323.45,421.99,314.69,421.21L314.69,421.21z"/>\n')
    out.write('<path class="st18" d="M439.29,333.39c4.48,18.78,6.56,37.15,6,54.19c-1.22,37.18-14.92,64.98-38.06,77.28\n')
    out.write(' c-23.14,12.31-53.85,8.11-85.35-11.66c-4.25-2.67-8.47-5.59-12.64-8.75l0,0c10.68-6.69,21.2-14.39,31.34-22.96\n')
    out.write(' c2.46-0.1,4.89-0.26,7.29-0.46c37.06-3.17,64.11-18.3,75.19-42.05c5.3-11.36,6.74-24.25,4.46-37.93c4.05-2.45,7.95-5,11.69-7.64\n')
    out.write(' L439.29,333.39z"/>\n')
    out.write('<path class="st19" d="M362.5,401.12c13.65-14.13,25.66-29.35,35.65-44.99c10.45-4.46,20.33-9.49,29.46-15.01\n')
    out.write(' c2.28,13.69,0.84,26.57-4.46,37.93c-11.08,23.76-38.13,38.88-75.19,42.05c-2.4,0.21-4.84,0.36-7.29,0.46\n')
    out.write(' C348.19,415.21,355.5,408.37,362.5,401.12z"/>\n')
    out.write('<path class="st20" d="M382.96,267.09c9.14,8.82,17.17,18.08,23.84,27.57c4.87,6.93,8.93,13.85,12.17,20.7\n')
    out.write(' c-5.33,13.68-12.44,27.5-21.06,41c-13.82,5.89-28.65,10.79-44.08,14.52c13.29-25.4,22.82-55.72,27.31-87.69\n')
    out.write(' C381.9,277.82,382.5,272.44,382.96,267.09L382.96,267.09z"/>\n')
    out.write('<path class="st21" d="M398.08,356.17c8.62-13.51,15.73-27.32,21.06-41c4.19,8.85,7,17.57,8.4,25.99\n')
    out.write(' C418.41,346.68,408.53,351.71,398.08,356.17z"/>\n')
    out.write('<path class="st22" d="M426.39,293.28c5.46,13.33,9.82,26.78,12.98,40.03l-0.1,0.02c-3.74,2.64-7.64,5.19-11.69,7.64\n')
    out.write(' c-1.4-8.42-4.21-17.14-8.4-25.99c2.83-7.26,5.15-14.47,6.94-21.59L426.39,293.28z"/>\n')
    out.write('<path class="st23" d="M425.79,293.57c-1.79,7.12-4.12,14.34-6.94,21.59c-3.24-6.85-7.3-13.78-12.17-20.7\n')
    out.write(' c-6.68-9.5-14.7-18.75-23.84-27.57l0,0c1.35-15.78,1.41-31.37,0.22-46.34c11.27,13.92,21.57,29.49,30.45,46.18\n')
    out.write(' c4.67,8.79,8.87,17.74,12.56,26.74L425.79,293.57z"/>\n')
    out.write('<path class="st24" d="M354.88,191.41c9.83,8.65,19.3,18.58,28.16,29.53c1.19,14.97,1.13,30.56-0.22,46.34l0,0\n')
    out.write(' c-20.15-19.44-45.7-36.8-74.13-50.06c-5.49-2.56-11.03-4.94-16.59-7.14l0,0c11.09-5.6,22.22-10.06,33.14-13.22\n')
    out.write(' c10.36-2.99,20.3-4.75,29.65-5.28L354.88,191.41z"/>\n')
    out.write('<path class="st25" d="M354.98,191.52c-9.35,0.54-19.29,2.29-29.65,5.28c-10.93,3.16-22.05,7.61-33.14,13.22l0,0\n')
    out.write(' c-15.7-6.19-31.6-10.91-47.19-14.05l0,0c17.54-6.97,36.54-12.35,56.24-15.83c11.42-2.01,22.84-3.34,34.1-4l0,0\n')
    out.write(' c6.66,4.48,13.22,9.57,19.62,15.2L354.98,191.52z"/>\n')
    out.write('<path class="st26" d="M426.21,293.66c-3.69-9-7.88-17.95-12.56-26.74c-8.88-16.69-19.18-32.26-30.45-46.18\n')
    out.write(' c-0.76-9.53-2.02-18.8-3.79-27.72l0,0c12.33,2.59,22.93,7.94,31.25,15.96c18.86,18.21,24.35,48.71,15.28,84.78L426.21,293.66z"/>\n')
    out.write('<path class="st27" d="M354.99,191.59c8.68-0.5,16.84,0.05,24.37,1.63l0,0c1.76,8.91,3.03,18.19,3.79,27.72\n')
    out.write(' c-8.87-10.95-18.33-20.88-28.16-29.53"/>\n')
    out.write('<path class="st28" d="M355.04,191.27c-6.41-5.63-12.97-10.72-19.62-15.2l0,0c13.7-0.8,27.16-0.6,40.08,0.59\n')
    out.write(' c1.51,5.31,2.82,10.8,3.93,16.42l0,0c-7.53-1.58-15.69-2.13-24.37-1.63L355.04,191.27z"/>\n')
    out.write('<path class="st29" d="M217.63,192.09c6.49-13.7,15.9-24.13,27.8-30.46c23.14-12.31,53.85-8.11,85.35,11.66\n')
    out.write(' c1.47,0.93,2.95,1.88,4.41,2.87l0,0c-11.26,0.66-22.68,1.99-34.1,4c-19.7,3.47-38.71,8.86-56.24,15.83l0,0\n')
    out.write(' C235.64,194.13,226.53,192.82,217.63,192.09L217.63,192.09z"/>\n')
    out.write('<path class="st30" d="M217.51,192.09c-5.32-0.44-10.56-0.67-15.71-0.69l0,0c7.73-19.66,17.5-37.36,28.86-51.96\n')
    out.write(' c22.84-29.36,50.26-43.8,76.22-40.15c25.96,3.65,48.33,25.09,62.2,59.61c2.28,5.69,4.31,11.65,6.07,17.85\n')
    out.write(' c-12.92-1.18-26.38-1.38-40.08-0.59l0,0c-1.47-0.99-2.94-1.94-4.41-2.87c-31.51-19.77-62.21-23.97-85.35-11.66\n')
    out.write(' C233.4,167.96,223.99,178.39,217.51,192.09L217.51,192.09z"/>\n')
    out.write('<path class="st31" d="M426.01,293.45c9.07-36.07,3.57-66.57-15.28-84.78c-8.31-8.03-18.92-13.38-31.25-15.96l0,0\n')
    out.write(' c-1.11-5.63-2.42-11.11-3.93-16.42c18.21,1.67,35.35,5.29,50.61,10.8c34.98,12.65,57.19,34.27,61.74,60.08s-8.93,53.72-37.47,77.57\n')
    out.write(' c-3.57,2.98-7.33,5.87-11.28,8.66l0.1-0.02c-3.16-13.25-7.51-26.7-12.98-40.03L426.01,293.45z"/>\n')

    # Gene names
    out.write(f'<text transform="matrix(1 0 0 1 66.0597 193.0202)" class="st32 st33">{gene_list[0]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 256.4101 90.4796)" class="st32 st33">{gene_list[1]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 436.3701 184.3349)" class="st32 st33">{gene_list[2]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 369.0695 498.2953)" class="st32 st33">{gene_list[3]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 161.0899 498.2954)" class="st32 st33">{gene_list[4]}</text>\n')

    # Counts
    out.write('<g>\n')
    count_text_str = '\n</g>\n<g>\n'.join(count_text_out)
    out.write(f'{count_text_str}\n')
    out.write('</g>\n')

    # Gradient
    out.write('<linearGradient id="SVGID_1_" gradientUnits="userSpaceOnUse" x1="37.5" y1="525" x2="37.5" y2="410.9666">\n')
    out.write('\t<stop  offset="0" style="stop-color:#FFFFFF"/>\n')
    out.write('\t<stop  offset="1" style="stop-color:#FF0000"/>\n')
    out.write('</linearGradient>\n')
    out.write('<rect x="27" y="411" class="st36" width="21" height="114"/>\n')

    # Rectangle around gradient
    out.write('<rect x="27" y="411" class="st37" width="21" height="114"/>\n')
    out.write('<line class="st1" x1="54.4" y1="411" x2="48" y2="411"/>\n')
    out.write('<line class="st1" x1="54.4" y1="439.5" x2="48" y2="439.5"/>\n')
    out.write('<line class="st1" x1="54.4" y1="468" x2="48" y2="468"/>\n')
    out.write('<line class="st1" x1="54.4" y1="496.5" x2="48" y2="496.5"/>\n')
    out.write('<line class="st1" x1="54.4" y1="525" x2="48" y2="525"/>\n')
    
    # Text for legend
    out.write('<g>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 415)" class="st34 st35">{"{:.1%}".format(max_per)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 443.5)" class="st34 st35">{"{:.1%}".format(max_per/4*3)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 472)" class="st34 st35">{"{:.1%}".format(max_per/4*2)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 500.5)" class="st34 st35">{"{:.1%}".format(max_per/4*1)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 529)" class="st34 st35">{"{:.1%}".format(0)}</text>\n')
    out.write('\t<text transform="matrix(6.123234e-17 1 -1 6.123234e-17 100 413.9671)" class="st34 st35">Percentage of peaks</text>\n')
    out.write('</g>\n')

    out.write('</svg>\n')

    out.close()

# Make the petal plot
def petal_three_col(infile_name, outfile_name, min_per, max_per, c1, c2, c3, gene_list):
    
    color_resolution = 1000
    max_factor = 1 / max_per  # e.g. 1/0.406
    min_factor = 1 / min_per  # e.g. 1/-0.406

    c1_rgb = rgb = matplotlib.colors.cnames[c1]
    c2_rgb = rgb = matplotlib.colors.cnames[c2]
    c3_rgb = rgb = matplotlib.colors.cnames[c3]

    out = open(outfile_name, 'w')
    table = open(infile_name, 'r')

    out.write('<?xml version="1.0" encoding="utf-8"?>\n')
    out.write('<!-- Generator: Adobe Illustrator 25.1.0, SVG Export Plug-In . SVG Version: 6.00 Build 0)  -->\n')
    out.write('<svg version="1.1" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px"\n')
    out.write('  viewBox="0 0 576 576" style="enable-background:new 0 0 576 576;" xml:space="preserve">\n')
    out.write('<style type="text/css">\n')
    out.write('\t.st1{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')

    # Initialize code dictionary
    combo = ['10000', '01000', '00100', '00010', '00001', '11000', '10100', '10010', '10001', '01100', '01010', '01001', '00110', 
    '00101', '00011', '11100', '11010', '11001', '10110', '10101', '10011', '01110', '01101', '01011', '00111', '11110', '11101',
    '11011', '10111', '01111', '11111']
    code_to_count = {}
    code_to_color = {}
    for code in combo:
        code_to_count[code] = 0
        code_to_color[code] = c2

    # Modify values from table
    count_text_out = []
    for line in table:
        tab = line.rstrip().split('\t')
        code = tab[1]
        code_to_count[code] = int(tab[4])
        
        color_out = c2
        color_input_value = float(tab[5])
        if color_input_value < 0:
            # Second color to first color
            color_frac = color_input_value * min_factor
            color_out = colorFader(c2,c1,color_frac)
            code_to_color[code] = color_out
        else:
            # Second color to third color
            color_frac = color_input_value * max_factor
            color_out = colorFader(c2,c3,color_frac)
            code_to_color[code] = color_out

    # Color and counts
    for code in combo:
        color_code = code_to_color[code]
        count = code_to_count[code]
        if code == '10000':
            out.write('\t.st2{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 120.5335 248.7162)" class="st34 st35">{count}</text>')
        elif code == '01000':
            out.write('\t.st30{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 278.6256 133.2766)" class="st34 st35">{count}</text>')
        elif code == '00100':
            out.write('\t.st31{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 441.5117 248.7159)" class="st34 st35">{count}</text>')
        elif code == '00010':
            out.write('\t.st18{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 369.7153 447.2975)" class="st34 st35">{count}</text>')
        elif code == '00001':
            out.write('\t.st15{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 187.6006 447.2975)" class="st34 st35">{count}</text>')
        elif code == '11000':
            out.write('\t.st6{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 186.6005 212.9799)" class="st34 st35">{count}</text>')
        elif code == '10100':
            out.write('\t.st3{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 148.0147 303.1999)" class="st34 st35">{count}</text>')
        elif code == '10010':
            out.write('\t.st19{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 381.3955 389.2256)" class="st34 st35">{count}</text>')
        elif code == '10001':
            out.write('\t.st10{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 161.6805 363.3396)" class="st34 st35">{count}</text>')
        elif code == '01100':
            out.write('\t.st28{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 344.9078 185.8387)" class="st34 st35">{count}</text>')
        elif code == '01010':
            out.write('\t.st29{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 246.6725 179.452)" class="st34 st35">{count}</text>')
        elif code == '01001':
            out.write('\t.st12{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 232.0956 424.2198)" class="st34 st35">{count}</text>')
        elif code == '00110':
            out.write('\t.st22{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 412.5907 323.2693)" class="st34 st35">{count}</text>')
        elif code == '00101':
            out.write('\t.st26{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 385.5107 225.0196)" class="st34 st35">{count}</text>')
        elif code == '00011':
            out.write('\t.st16{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 294.0507 438.8972)" class="st34 st35">{count}</text>')
        elif code == '11100':
            out.write('\t.st4{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 179.8049 259.4544)" class="st34 st35">{count}</text>')
        elif code == '11010':
            out.write('\t.st5{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 211.1205 202.4147)" class="st34 st35">{count}</text>')
        elif code == '11001':
            out.write('\t.st11{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 211.1205 386.3497)" class="st34 st35">{count}</text>')
        elif code == '10110':
            out.write('\t.st21{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 399.7697 342.9097)" class="st34 st35">{count}</text>')
        elif code == '10101':
            out.write('\t.st8{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 160.6805 342.9097)" class="st34 st35">{count}</text>')
        elif code == '10011':
            out.write('\t.st17{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 339.5356 396.0097)" class="st34 st35">{count}</text>')
        elif code == '01110':
            out.write('\t.st25{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 284.6205 195.9897)" class="st34 st35">{count}</text>')
        elif code == '01101':
            out.write('\t.st27{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 355.9945 204.5796)" class="st34 st35">{count}</text>')
        elif code == '01011':
            out.write('\t.st14{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 280.6556 427.3506)" class="st34 st35">{count}</text>')
        elif code == '00111':
            out.write('\t.st23{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 381.3955 267.8547)" class="st34 st35">{count}</text>')
        elif code == '11110':
            out.write('\t.st0{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 221.8954 226.0797)" class="st34 st35">{count}</text>')
        elif code == '11101':
            out.write('\t.st9{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 187.6006 331.6697)" class="st34 st35">{count}</text>')
        elif code == '11011':
            out.write('\t.st13{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 284.6205 402.5913)" class="st34 st35">{count}</text>')
        elif code == '10111':
            out.write('\t.st20{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 377.3163 325.96)" class="st34 st35">{count}</text>')
        elif code == '01111':
            out.write('\t.st24{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 331.438 220.9398)" class="st34 st35">{count}</text>')
        elif code == '11111':
            out.write('\t.st7{fill:'+str(color_code)+';stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
            count_text_out.append(f'\t<text transform="matrix(1 0 0 1 278.998 303.2002)" class="st34 st35">{count}</text>')
        

    out.write('\t.st32{font-family:\'Helvetica-Oblique\';}\n')
    out.write('\t.st33{font-size:24px;}\n')
    out.write('\t.st34{font-family:\'Helvetica\';}\n')
    out.write('\t.st35{font-size:12px;}\n')
    out.write('\t.st36{fill:none;stroke:#FFFFFF;stroke-width:0;stroke-miterlimit:10;}\n')
    out.write('\t.st37{fill:none;stroke:#FFFFFF;stroke-width:0;stroke-miterlimit:10;}\n')
    out.write('\t.st38{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}\n')
    out.write('</style>\n')

    # Locations
    out.write('<path class="st0" d="M209.86,275.16c-1.98-12.51-2.81-24.69-2.43-36.25c0.31-9.35,1.41-18.11,3.24-26.18l0,0\n')
    out.write(' c10.53-6.27,22.03-11.89,34.24-16.74l0,0c15.59,3.14,31.49,7.86,47.19,14.05l0,0c-25.18,12.73-50.16,31.37-72.16,54.15\n')
    out.write(' C216.48,267.78,213.12,271.44,209.86,275.16L209.86,275.16z"/>\n')
    out.write('<rect class="st1" width="576" height="576"/>\n')
    out.write('<path class="st2" d="M161.04,355.52c-13.44-11.52-24.99-23.96-34.06-36.85c-21.4-30.43-27.19-60.87-16.12-84.63\n')
    out.write(' c11.08-23.76,38.13-38.88,75.19-42.05c5.12-0.44,10.38-0.64,15.74-0.62l0,0c-4.33,11.02-8.03,22.66-11,34.71l0,0\n')
    out.write(' c-1.56,1.21-3.08,2.44-4.56,3.68c-28.54,23.85-42.02,51.76-37.47,77.57c2.18,12.35,8.39,23.73,18.12,33.64l0,0\n')
    out.write(' C164.72,345.83,162.77,350.69,161.04,355.52L161.04,355.52z"/>\n')
    out.write('<path class="st3" d="M190.94,226.02c-2.39,9.68-4.31,19.62-5.73,29.73c-2.66,18.9-3.46,37.72-2.48,55.76l0,0\n')
    out.write(' c-6.05,9.71-11.32,19.56-15.7,29.4l0,0c-9.72-9.9-15.94-21.29-18.12-33.64c-4.55-25.81,8.93-53.72,37.47-77.57\n')
    out.write(' C187.86,228.46,189.38,227.23,190.94,226.02L190.94,226.02z"/>\n')
    out.write('<path class="st4" d="M210.64,212.98c-1.84,8.07-2.94,16.83-3.24,26.19c-0.38,11.57,0.45,23.75,2.43,36.25l0,0\n')
    out.write(' c-10.2,11.67-19.29,23.93-27.09,36.43l0,0c-0.98-18.04-0.17-36.86,2.48-55.76c1.42-10.11,3.34-20.05,5.73-29.73l0,0\n')
    out.write(' C197.03,221.63,203.62,217.16,210.64,212.98L210.64,212.98z"/>\n')
    out.write('<path class="st5" d="M210.61,212.73c1.71-7.5,4.06-14.41,7.01-20.63l0,0c8.94,0.73,18.09,2.05,27.36,3.92l0.02-0.08\n')
    out.write(' C232.73,200.79,221.18,206.43,210.61,212.73L210.61,212.73z"/>\n')
    out.write('<path class="st6" d="M210.62,212.7c-7.02,4.18-13.6,8.65-19.68,13.38l0,0c2.97-12.05,6.67-23.69,11-34.71l0,0\n')
    out.write(' c5.14,0.02,10.38,0.25,15.69,0.69l0,0C214.68,198.29,212.33,205.19,210.62,212.7L210.62,212.7z"/>\n')
    out.write('<path class="st7" d="M382.83,267.21c-0.46,5.35-1.06,10.73-1.82,16.11c-4.49,31.97-14.02,62.29-27.31,87.69l0,0\n')
    out.write(' c-5.99,1.45-12.07,2.72-18.22,3.8c-29.93,5.28-59.85,5.86-86.84,1.91l0,0c-3.39-5.47-6.62-11.11-9.69-16.87\n')
    out.write(' c-14.75-27.74-24.72-57.02-29.1-84.71l0,0c3.25-3.72,6.62-7.39,10.08-10.98c22-22.78,46.98-41.42,72.16-54.15l0,0\n')
    out.write(' c5.56,2.2,11.1,4.58,16.6,7.14C337.13,230.42,362.69,247.77,382.83,267.21L382.83,267.21z"/>\n')
    out.write('<path class="st8" d="M167.04,340.98c4.38-9.84,9.65-19.69,15.7-29.4l0,0c0.88,16.27,3.22,31.91,6.95,46.36l0,0\n')
    out.write(' C180.82,352.97,173.22,347.28,167.04,340.98L167.04,340.98z"/>\n')
    out.write('<path class="st9" d="M209.99,275.17c4.38,27.69,14.35,56.97,29.1,84.71c3.07,5.76,6.3,11.4,9.69,16.87l0,0\n')
    out.write(' c-13.51-1.97-26.28-5.08-37.95-9.3c-7.61-2.75-14.62-5.93-20.97-9.48l0,0c-3.73-14.46-6.07-30.09-6.95-46.36l0,0\n')
    out.write(' C190.7,299.1,199.79,286.84,209.99,275.17L209.99,275.17z"/>\n')
    out.write('<path class="st10" d="M166.93,340.91c-2.17,4.86-4.12,9.72-5.84,14.56c10.96,9.39,23.17,18.17,36.32,26.09\n')
    out.write(' c-0.19-0.45-0.37-0.91-0.56-1.36c-2.83-7.04-5.26-14.51-7.28-22.33l0,0C180.71,352.9,173.11,347.21,166.93,340.91L166.93,340.91z"/>\n')
    out.write('<path class="st11" d="M248.61,376.73c8.52,13.77,18,26.56,28.13,38.02c-17.22-4.33-34.65-10.52-51.62-18.43\n')
    out.write(' c-9.56-4.46-18.79-9.37-27.6-14.68c-0.19-0.45-0.37-0.91-0.56-1.36c-2.83-7.04-5.26-14.51-7.28-22.33l0,0\n')
    out.write(' c6.35,3.55,13.36,6.73,20.97,9.48C222.33,371.65,235.11,374.75,248.61,376.73L248.61,376.73z"/>\n')
    out.write('<path class="st12" d="M259.06,439.81c-25.61-3.6-47.74-24.53-61.64-58.24c8.81,5.31,18.04,10.22,27.6,14.68\n')
    out.write(' c16.97,7.91,34.4,14.11,51.62,18.43c6.02,6.81,12.27,13.15,18.69,18.94C283.32,439.32,271.03,441.49,259.06,439.81z"/>\n')
    out.write('<path class="st13" d="M314.85,421.14c-12.34-1.1-25.08-3.29-37.94-6.53c-10.13-11.46-19.61-24.24-28.13-38.02l0,0\n')
    out.write(' c26.98,3.94,56.9,3.36,86.83-1.91c6.15-1.08,12.23-2.35,18.22-3.8c-5.48,10.48-11.61,20.13-18.28,28.71\n')
    out.write(' C328.98,408.03,322.03,415.24,314.85,421.14L314.85,421.14z"/>\n')
    out.write('<path class="st14" d="M314.68,421.28c-6.27,5.15-12.72,9.31-19.26,12.41c-6.42-5.79-12.67-12.13-18.69-18.94\n')
    out.write(' C289.6,417.98,302.34,420.17,314.68,421.28L314.68,421.28z"/>\n')
    out.write('<path class="st15" d="M308.75,444.71c-17.25,10.81-34.94,18.99-52.13,23.96c-35.74,10.32-66.41,5.89-85.26-12.31\n')
    out.write(' c-18.86-18.21-24.35-48.71-15.28-84.78c1.35-5.36,3-10.77,4.94-16.22c10.96,9.39,23.17,18.17,36.32,26.09\n')
    out.write(' c13.9,33.72,36.03,54.64,61.64,58.24c11.97,1.68,24.26-0.48,36.26-6.19C299.67,437.51,304.18,441.25,308.75,444.71L308.75,444.71z"\n')
    out.write(' />\n')
    out.write('<path class="st16" d="M314.57,421.14c8.76,0.78,17.33,1.01,25.59,0.66c-10.13,8.57-20.65,16.27-31.34,22.96l0,0\n')
    out.write(' c-4.57-3.46-9.08-7.2-13.51-11.2C301.85,430.45,308.3,426.3,314.57,421.14L314.57,421.14z"/>\n')
    out.write('<path class="st17" d="M314.69,421.21c7.18-5.9,14.13-13.11,20.7-21.55c6.67-8.58,12.8-18.22,18.28-28.71\n')
    out.write(' c15.44-3.73,30.26-8.63,44.08-14.52c-9.98,15.64-22,30.86-35.65,44.99c-7,7.25-14.31,14.09-21.83,20.44\n')
    out.write(' C332.01,422.21,323.45,421.99,314.69,421.21L314.69,421.21z"/>\n')
    out.write('<path class="st18" d="M439.29,333.39c4.48,18.78,6.56,37.15,6,54.19c-1.22,37.18-14.92,64.98-38.06,77.28\n')
    out.write(' c-23.14,12.31-53.85,8.11-85.35-11.66c-4.25-2.67-8.47-5.59-12.64-8.75l0,0c10.68-6.69,21.2-14.39,31.34-22.96\n')
    out.write(' c2.46-0.1,4.89-0.26,7.29-0.46c37.06-3.17,64.11-18.3,75.19-42.05c5.3-11.36,6.74-24.25,4.46-37.93c4.05-2.45,7.95-5,11.69-7.64\n')
    out.write(' L439.29,333.39z"/>\n')
    out.write('<path class="st19" d="M362.5,401.12c13.65-14.13,25.66-29.35,35.65-44.99c10.45-4.46,20.33-9.49,29.46-15.01\n')
    out.write(' c2.28,13.69,0.84,26.57-4.46,37.93c-11.08,23.76-38.13,38.88-75.19,42.05c-2.4,0.21-4.84,0.36-7.29,0.46\n')
    out.write(' C348.19,415.21,355.5,408.37,362.5,401.12z"/>\n')
    out.write('<path class="st20" d="M382.96,267.09c9.14,8.82,17.17,18.08,23.84,27.57c4.87,6.93,8.93,13.85,12.17,20.7\n')
    out.write(' c-5.33,13.68-12.44,27.5-21.06,41c-13.82,5.89-28.65,10.79-44.08,14.52c13.29-25.4,22.82-55.72,27.31-87.69\n')
    out.write(' C381.9,277.82,382.5,272.44,382.96,267.09L382.96,267.09z"/>\n')
    out.write('<path class="st21" d="M398.08,356.17c8.62-13.51,15.73-27.32,21.06-41c4.19,8.85,7,17.57,8.4,25.99\n')
    out.write(' C418.41,346.68,408.53,351.71,398.08,356.17z"/>\n')
    out.write('<path class="st22" d="M426.39,293.28c5.46,13.33,9.82,26.78,12.98,40.03l-0.1,0.02c-3.74,2.64-7.64,5.19-11.69,7.64\n')
    out.write(' c-1.4-8.42-4.21-17.14-8.4-25.99c2.83-7.26,5.15-14.47,6.94-21.59L426.39,293.28z"/>\n')
    out.write('<path class="st23" d="M425.79,293.57c-1.79,7.12-4.12,14.34-6.94,21.59c-3.24-6.85-7.3-13.78-12.17-20.7\n')
    out.write(' c-6.68-9.5-14.7-18.75-23.84-27.57l0,0c1.35-15.78,1.41-31.37,0.22-46.34c11.27,13.92,21.57,29.49,30.45,46.18\n')
    out.write(' c4.67,8.79,8.87,17.74,12.56,26.74L425.79,293.57z"/>\n')
    out.write('<path class="st24" d="M354.88,191.41c9.83,8.65,19.3,18.58,28.16,29.53c1.19,14.97,1.13,30.56-0.22,46.34l0,0\n')
    out.write(' c-20.15-19.44-45.7-36.8-74.13-50.06c-5.49-2.56-11.03-4.94-16.59-7.14l0,0c11.09-5.6,22.22-10.06,33.14-13.22\n')
    out.write(' c10.36-2.99,20.3-4.75,29.65-5.28L354.88,191.41z"/>\n')
    out.write('<path class="st25" d="M354.98,191.52c-9.35,0.54-19.29,2.29-29.65,5.28c-10.93,3.16-22.05,7.61-33.14,13.22l0,0\n')
    out.write(' c-15.7-6.19-31.6-10.91-47.19-14.05l0,0c17.54-6.97,36.54-12.35,56.24-15.83c11.42-2.01,22.84-3.34,34.1-4l0,0\n')
    out.write(' c6.66,4.48,13.22,9.57,19.62,15.2L354.98,191.52z"/>\n')
    out.write('<path class="st26" d="M426.21,293.66c-3.69-9-7.88-17.95-12.56-26.74c-8.88-16.69-19.18-32.26-30.45-46.18\n')
    out.write(' c-0.76-9.53-2.02-18.8-3.79-27.72l0,0c12.33,2.59,22.93,7.94,31.25,15.96c18.86,18.21,24.35,48.71,15.28,84.78L426.21,293.66z"/>\n')
    out.write('<path class="st27" d="M354.99,191.59c8.68-0.5,16.84,0.05,24.37,1.63l0,0c1.76,8.91,3.03,18.19,3.79,27.72\n')
    out.write(' c-8.87-10.95-18.33-20.88-28.16-29.53"/>\n')
    out.write('<path class="st28" d="M355.04,191.27c-6.41-5.63-12.97-10.72-19.62-15.2l0,0c13.7-0.8,27.16-0.6,40.08,0.59\n')
    out.write(' c1.51,5.31,2.82,10.8,3.93,16.42l0,0c-7.53-1.58-15.69-2.13-24.37-1.63L355.04,191.27z"/>\n')
    out.write('<path class="st29" d="M217.63,192.09c6.49-13.7,15.9-24.13,27.8-30.46c23.14-12.31,53.85-8.11,85.35,11.66\n')
    out.write(' c1.47,0.93,2.95,1.88,4.41,2.87l0,0c-11.26,0.66-22.68,1.99-34.1,4c-19.7,3.47-38.71,8.86-56.24,15.83l0,0\n')
    out.write(' C235.64,194.13,226.53,192.82,217.63,192.09L217.63,192.09z"/>\n')
    out.write('<path class="st30" d="M217.51,192.09c-5.32-0.44-10.56-0.67-15.71-0.69l0,0c7.73-19.66,17.5-37.36,28.86-51.96\n')
    out.write(' c22.84-29.36,50.26-43.8,76.22-40.15c25.96,3.65,48.33,25.09,62.2,59.61c2.28,5.69,4.31,11.65,6.07,17.85\n')
    out.write(' c-12.92-1.18-26.38-1.38-40.08-0.59l0,0c-1.47-0.99-2.94-1.94-4.41-2.87c-31.51-19.77-62.21-23.97-85.35-11.66\n')
    out.write(' C233.4,167.96,223.99,178.39,217.51,192.09L217.51,192.09z"/>\n')
    out.write('<path class="st31" d="M426.01,293.45c9.07-36.07,3.57-66.57-15.28-84.78c-8.31-8.03-18.92-13.38-31.25-15.96l0,0\n')
    out.write(' c-1.11-5.63-2.42-11.11-3.93-16.42c18.21,1.67,35.35,5.29,50.61,10.8c34.98,12.65,57.19,34.27,61.74,60.08s-8.93,53.72-37.47,77.57\n')
    out.write(' c-3.57,2.98-7.33,5.87-11.28,8.66l0.1-0.02c-3.16-13.25-7.51-26.7-12.98-40.03L426.01,293.45z"/>\n')

    # Gene names
    out.write(f'<text transform="matrix(1 0 0 1 66.0597 193.0202)" class="st32 st33">{gene_list[0]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 256.4101 90.4796)" class="st32 st33">{gene_list[1]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 436.3701 184.3349)" class="st32 st33">{gene_list[2]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 369.0695 498.2953)" class="st32 st33">{gene_list[3]}</text>\n')
    out.write(f'<text transform="matrix(1 0 0 1 161.0899 498.2954)" class="st32 st33">{gene_list[4]}</text>\n')

    # Counts
    out.write('<g>\n')
    count_text_str = '\n</g>\n<g>\n'.join(count_text_out)
    out.write(f'{count_text_str}\n')
    out.write('</g>\n')

    # Gradient part 1
    out.write('<linearGradient id="SVGID_1_" gradientUnits="userSpaceOnUse" x1="37.5" y1="468" x2="37.5" y2="411">\n')
    out.write(f'\t<stop  offset="0" style="stop-color:{c2_rgb}"/>\n')
    out.write(f'\t<stop  offset="1" style="stop-color:{c3_rgb}"/>\n')
    out.write('</linearGradient>\n')
    out.write('<rect x="27" y="411" class="st36" width="21" height="57"/>\n')

    # Gradient part 2
    out.write('<linearGradient id="SVGID_2_" gradientUnits="userSpaceOnUse" x1="37.5" y1="525" x2="37.5" y2="468">\n')
    out.write(f'\t<stop  offset="0" style="stop-color:{c1_rgb}"/>\n')
    out.write(f'\t<stop  offset="1" style="stop-color:{c2_rgb}"/>\n')
    out.write('</linearGradient>\n')
    out.write('<rect x="27" y="468" class="st37" width="21" height="57"/>\n')

    # Rectangle around gradient
    out.write('<rect x="27" y="411" class="st38" width="21" height="114"/>\n')
    out.write('<line class="st1" x1="54.4" y1="411" x2="48" y2="411"/>\n')
    out.write('<line class="st1" x1="54.4" y1="439.5" x2="48" y2="439.5"/>\n')
    out.write('<line class="st1" x1="54.4" y1="468" x2="48" y2="468"/>\n')
    out.write('<line class="st1" x1="54.4" y1="496.5" x2="48" y2="496.5"/>\n')
    out.write('<line class="st1" x1="54.4" y1="525" x2="48" y2="525"/>\n')
    
    # Text for legend
    out.write('<g>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 415)" class="st34 st35">{"{:.2f}".format(max_per)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 443.5)" class="st34 st35">{"{:.2f}".format(max_per/2)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.6576 472)" class="st34 st35">{"{:.0f}".format(0)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.2856 500.5)" class="st34 st35">{"{:.2f}".format(min_per/2)}</text>\n')
    out.write(f'\t<text transform="matrix(1 0 0 1 57.2856 529)" class="st34 st35">{"{:.2f}".format(min_per)}</text>\n')
    out.write('\t<text transform="matrix(6.123234e-17 1 -1 6.123234e-17 100 438)" class="st34 st35">Cohen\'s h</text>\n')
    out.write('</g>\n')

    out.write('</svg>\n')

    out.close()


if __name__ == '__main__':
    main(
        *parse_args()
    )