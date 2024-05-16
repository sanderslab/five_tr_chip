# Code for assessing overlap between replicates
# Example usage:
# python compare_replicates.py -a arid1b_peaks_hg38.txt -b atac_peaks_hg38.txt -mo 0.01

# Requires bedtools to be installed

import csv
import subprocess
from argparse import ArgumentParser
from matplotlib import pyplot as plt

def __key_input_bed(in_file):
    """
    Add a column to the bed file to reference the region's
    original line number
    """
    line_count = 1
    anno_file = f'{in_file}.linekey.bed'
    with open(anno_file, 'w') as out_bed:
        w = csv.writer(out_bed, delimiter='\t')
        with open(in_file, 'r') as in_bed:
            for line in in_bed:
                fields = [i.strip() for i in line.split('\t')]
                fields.append(line_count)
                w.writerow(fields)
                line_count += 1
    return anno_file, line_count


def __run_bedtools_intersect(in_file_a, in_file_b, label_a, label_b, min_overlap, working_dir):
    """
    Bedtools intersect output cols
    A_chrom A_start A_end = fields[:3]
    A_og_TF_line = fields[10]
    B_chrom B_start B_end = fields[11:14]
    B_og_TF_line = fields[21]
    num_bases_overlap = fields[22]
    """
    cmd = ['bedtools', 'intersect', '-a', in_file_a, '-b', in_file_b, '-wao', '-r', '-f', str(min_overlap)]
    outfile = f'{working_dir}/{label_a}_{label_b}_{min_overlap}_intersect.bed'
    print(cmd)
    with open(outfile, 'w') as fh:
        subprocess.call(cmd, stdout=fh)
    return outfile


def __make_overlap_dict(intersect_bed, key_by):
    """
    Make a dict from the bedtools intersect file key'd by
    A's line numbers.

    overlap_dict[key]['liftover_region'] is the region lifted over, even
    if that region isnt a TF in the other species; take this from
    the raw lifted over file.
    """
    overlap_dict = dict()
    with open(intersect_bed, 'r') as in_bed:
        for line in in_bed:
            fields = [i.strip() for i in line.split('\t')]
            # print(fields)
            a_chrom, a_start, a_end = fields[:3]
            a_og_tf_line = fields[3]
            b_chrom, b_start, b_end = fields[4:7]
            b_og_tf_line = fields[7]
            num_bases_overlap = fields[8]
            # a_og_tf_line = fields[10]
            # b_chrom, b_start, b_end = fields[11:14]
            # b_og_tf_line = fields[21]
            # num_bases_overlap = fields[22]

            # Only add regions that had an intersect to overlap_dict
            if num_bases_overlap != '0':
                if key_by == 'A':
                    if a_og_tf_line not in overlap_dict:
                        overlap_dict[a_og_tf_line] = list()
                        overlap_dict[a_og_tf_line].append((b_chrom, b_start, b_end, num_bases_overlap))
                    else:
                        # raise ValueError(f'Multiple TF overlaps for A line {A_og_TF_line}?')
                        overlap_dict[a_og_tf_line].append((b_chrom, b_start, b_end, num_bases_overlap))
                elif key_by == 'B':
                    if b_og_tf_line not in overlap_dict:
                        overlap_dict[b_og_tf_line] = list()
                        overlap_dict[b_og_tf_line].append((a_chrom, a_start, a_end, num_bases_overlap))
                    else:
                        # raise ValueError(f'Multiple TF overlaps for A line {A_og_TF_line}?')
                        overlap_dict[b_og_tf_line].append((a_chrom, a_start, a_end, num_bases_overlap))

    return overlap_dict


def __anno_overlaps(infile, overlap_dict, a_file, b_file, a_prefix, b_prefix, min_overlap, working_dir):
    """
    Add columns to the original TF bed file for:
        - lifted over region if liftover was successful
        - TF region from the other species if lifted over region is a TF
    Use the line-key'd bed file as input
    """
    # overlapping_TFs = 0
    anno_file = f'{working_dir}/{a_prefix}_{b_prefix}.overlap_anno.{min_overlap}.bed'
    print(f'Writing to: {anno_file}')
    with open(anno_file, 'w') as out_bed:
        w = csv.writer(out_bed, delimiter='\t')
        header1 = [f'#Input1:{a_file}', f'Input2:{b_file}', f'min_overlap:{min_overlap}']
        header2 = [
            '#A_chrom', 'A_start', 'A_end',
            'B_overlaps (chrom:start-end)',
            '%ReciprocalOverlap', '%OverlapA', '%OverlapB', 'NumBasesOverlap']
        w.writerow(header1)
        w.writerow(header2)

        with open(infile, 'r') as in_bed:
            for line in in_bed:
                fields = [i.strip() for i in line.split('\t')]
                a_chrom, a_start, a_end = fields[:3]
                key = fields[-1]
                if key not in overlap_dict:
                    b_chrom, b_start, b_end, overlaps = ('.', '.', '.', '.')
                    recip_percent_overlap, a_percent_overlap, b_percent_overlap = ('.', '.', '.')
                    num_bases_overlap = '.'
                else:
                    overlap_list = overlap_dict[key]
                    overlaps = list()
                    recip_p_overlaps = list()
                    a_p_overlaps = list()
                    b_p_overlaps = list()
                    n_bp_overlaps = list()
                    for region in overlap_list:
                        b_chrom, b_start, b_end = region[:3]
                        num_bases_overlap = region[3]

                        a_percentage_overlap = str(
                            round((int(num_bases_overlap) / (int(a_end) - int(a_start))) * 100, 1))
                        b_percentage_overlap = str(
                            round((int(num_bases_overlap) / (int(b_end) - int(b_start))) * 100, 1))

                        recip_p_overlap = str(min(float(a_percentage_overlap), float(b_percentage_overlap)))

                        overlaps.append(f'{b_chrom}:{b_start}-{b_end}')
                        recip_p_overlaps.append(recip_p_overlap)
                        a_p_overlaps.append(a_percentage_overlap)
                        b_p_overlaps.append(b_percentage_overlap)
                        n_bp_overlaps.append(num_bases_overlap)

                    overlaps = ','.join(overlaps)
                    num_bases_overlap = ','.join(n_bp_overlaps)
                    recip_percent_overlap = ','.join(recip_p_overlaps)
                    a_percent_overlap = ','.join(a_p_overlaps)
                    b_percent_overlap = ','.join(b_p_overlaps)

                anno_row = fields[:-1] + [overlaps,
                                          recip_percent_overlap, a_percent_overlap, b_percent_overlap,
                                          num_bases_overlap]
                w.writerow(anno_row)

    return anno_file


def __plot_overlap_dist_hist(overlap_bed_a, overlap_bed_b, a_prefix, b_prefix, min_overlap, working_dir):
    """
    Plot a histogram of the %overlap of the two species' TFs (7th col)
    """
    overlaps_a = list()
    overlaps_b = list()
    total_overlap_bases_a = total_bases_a = total_nonoverlap_bases_a = 0
    total_peaks_a = overlap_peaks_a = nonoverlap_peaks_a = 0
    total_overlap_bases_b = total_bases_b = total_nonoverlap_bases_b = 0
    total_peaks_b = overlap_peaks_b = nonoverlap_peaks_b = 0

    with open(overlap_bed_a, 'r') as in_f:
        for line in in_f:
            if not line.startswith('#'):
                total_peaks_a += 1
                fields = [i.strip() for i in line.split('\t')]
                region_size_a = int(fields[2]) - int(fields[1])
                total_bases_a += region_size_a
                recip_p_overlap = fields[4]
                bases_overlap = fields[7]
                non_overlap_a = region_size_a
                if recip_p_overlap != '.':
                    # Overlaps at least one region in file B
                    overlap_peaks_a += 1
                    recip_p_overlap_list = recip_p_overlap.split(',')
                    overlaps_a.extend([float(i) for i in recip_p_overlap_list])
                    bases_overlap_list = bases_overlap.split(',')
                    for base_count in bases_overlap_list:
                        total_overlap_bases_a += int(base_count)
                        non_overlap_a -= int(base_count)

                    total_nonoverlap_bases_a += non_overlap_a
                else:
                    nonoverlap_peaks_a += 1
                    total_nonoverlap_bases_a += region_size_a

    plt.figure(figsize=(10, 6))
    plt.rcParams['pdf.use14corefonts'] = True
    plt.hist(overlaps_a, bins=101)
    plt.ylabel('Counts')
    plt.xlabel('Percent overlap')
    plt.title(f'Percent overlap distribution from File A: {a_prefix}\nto File B: {b_prefix}\n'
              f'w/ min_overlap={min_overlap}, total/overlap_peaks_A={total_peaks_a}/{overlap_peaks_a}, '
              f'total/overlap_peaks_B={total_peaks_b}/{overlap_peaks_b}')
    plt.savefig(f'{working_dir}/{a_prefix}_to_{b_prefix}_recip_percentage_overlap.hist.pdf')

    with open(overlap_bed_b, 'r') as in_f:
        for line in in_f:
            if not line.startswith('#'):
                total_peaks_b += 1
                fields = [i.strip() for i in line.split('\t')]
                region_size_b = int(fields[2]) - int(fields[1])
                total_bases_b += region_size_b
                recip_p_overlap = fields[4]
                bases_overlap = fields[7]
                if recip_p_overlap != '.':
                    # Overlaps at least one region in file A
                    overlap_peaks_b += 1
                    recip_p_overlap_list = recip_p_overlap.split(',')
                    overlaps_b.extend([float(i) for i in recip_p_overlap_list])
                    bases_overlap_list = bases_overlap.split(',')
                    non_overlap_b = region_size_b
                    for base_count in bases_overlap_list:
                        total_overlap_bases_b += int(base_count)
                        non_overlap_b -= int(base_count)

                    total_nonoverlap_bases_b += non_overlap_b
                else:
                    nonoverlap_peaks_b += 1
                    total_nonoverlap_bases_b += region_size_b

    a_per_peaks = f'{round(overlap_peaks_a / total_peaks_a * 100, 2)}%'
    a_per_bases = f'{round(total_overlap_bases_a / total_bases_a * 100, 2)}%'
    b_per_peaks = f'{round(overlap_peaks_b / total_peaks_b * 100, 2)}%'
    b_per_bases = f'{round(total_overlap_bases_b / total_bases_b * 100, 2)}%'

    plt.figure(figsize=(10, 6))
    plt.rcParams['pdf.use14corefonts'] = True
    plt.hist(overlaps_b, bins=101)
    plt.ylabel('Counts')
    plt.xlabel('Percent overlap')
    plt.title(f'Percent overlap distribution from File B: {b_prefix}\nto File A: {a_prefix}\n'
              f'w/ min_overlap={min_overlap}, total/overlap_peaks_B={total_peaks_b}/{overlap_peaks_b}, '
              f'total/overlap_peaks_A={total_peaks_a}/{overlap_peaks_a}')
    plt.savefig(f'{working_dir}/{b_prefix}_to_{a_prefix}_recip_percentage_overlap.hist.pdf')

    out_file_name = f'{working_dir}/{a_prefix}_to_{b_prefix}_summary.txt'
    print(f'Writing summary to {out_file_name}')
    out = open(out_file_name, 'w')
    out.write('File_A\tTotal_peaks_A\tOverlap_peaks_A\tNon-overlap_peaks_A\t'
              'Total_bases_A\tOverlap_bases_A\tNon-overlap_bases_A\t'
              'A_%_peaks_overlap\tA_%_bases_overlap\t'
              'File_B\tTotal_peaks_B\tOverlap_peaks_B\tNon-overlap_peaks_B\t'
              'Total_bases_B\tOverlap_bases_B\tNon-overlap_bases_B\t'
              'B_%_peaks_overlap\tB_%_bases_overlap\n')
    out.write(f'{a_prefix}\t{total_peaks_a}\t{overlap_peaks_a}\t{nonoverlap_peaks_a}\t'
              f'{total_bases_a}\t{total_overlap_bases_a}\t{total_nonoverlap_bases_a}\t'
              f'{a_per_peaks}\t{a_per_bases}\t'
              f'{b_prefix}\t{total_peaks_b}\t{overlap_peaks_b}\t{nonoverlap_peaks_b}\t'
              f'{total_bases_b}\t{total_overlap_bases_b}\t{total_nonoverlap_bases_b}\t'
              f'{b_per_peaks}\t{b_per_bases}\n')


def main(file_a, file_b, min_overlap):
    # curr_dir_contents = os.listdir('.')

    file_a_simp = file_a.split('/')[-1]
    file_b_simp = file_b.split('/')[-1]
    working_dir = '/'.join(file_a.split("/")[:-1])
    a_prefix = file_a_simp.split('.')[0]
    b_prefix = file_b_simp.split('.')[0]

    # Annotate line numbers to original beds
    print('Adding line number keys to bed')
    a_key_file, a_num_regions = __key_input_bed(file_a)
    b_key_file, b_num_regions = __key_input_bed(file_b)

    # Intersect lifted over/merged human TFs (now in mm coords) w/ original mm TFs
    print('Running bedtools intersect')
    intersect_bed = __run_bedtools_intersect(a_key_file, b_key_file, a_prefix, b_prefix, min_overlap, working_dir)

    # Make overlap dicts for annotation
    a_overlap_dict = __make_overlap_dict(
        intersect_bed=intersect_bed, key_by='A')

    b_overlap_dict = __make_overlap_dict(
        intersect_bed=intersect_bed, key_by='B')

    # Annotate original hg and mm TF bed files with overlap information
    a_anno_overlap = __anno_overlaps(
        infile=a_key_file,
        overlap_dict=a_overlap_dict,
        a_file=a_key_file,
        b_file=b_key_file,
        a_prefix=a_prefix,
        b_prefix=b_prefix,
        min_overlap=min_overlap,
        working_dir=working_dir)

    b_anno_overlap = __anno_overlaps(
        infile=b_key_file,
        overlap_dict=b_overlap_dict,
        a_file=b_key_file,
        b_file=a_key_file,
        a_prefix=b_prefix,
        b_prefix=a_prefix,
        min_overlap=min_overlap,
        working_dir=working_dir)

    __plot_overlap_dist_hist(
        overlap_bed_a=a_anno_overlap,
        overlap_bed_b=b_anno_overlap,
        a_prefix=a_prefix,
        b_prefix=b_prefix,
        min_overlap=min_overlap,
        working_dir=working_dir)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-a', '--file-a', type=str, required=True)
    parser.add_argument('-b', '--file-b', type=str, required=True)
    parser.add_argument('-mo', '--min-overlap', type=str, required=True)
    args = parser.parse_args()
    return args.file_a, args.file_b, float(args.min_overlap)


if __name__ == '__main__':
    main(
        *parse_args()
    )
