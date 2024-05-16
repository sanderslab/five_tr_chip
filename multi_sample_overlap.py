# Function for the assessment of overlap between multiple files of ChIP-seq data (deployed in Jupyter Lab)
# Example usage:
# spearman_peaks(prom_liver_brain_list, 10000, 'prom_liver_brain_cor_10000.txt')
# prom_liver_brain_list is a list of files of ChIP-seq peaks sorted by log10P 
# 10000 is the number of peaks to use across the peak files
# prom_liver_brain_cor_10000.txt is the output file name

import os.path
from os import path
from scipy import stats
from scipy.stats import uniform, randint
from scipy.stats import spearmanr  # Soearman test

# Liver overlap, function to assess correlation between files, Fig 2B/C

def spearman_peaks(file_list, max_peak_count, out_file):
    out = open(out_file, 'w')
    out.write('File_A_Path\tFile_B_Path\tFile_A\tFile_B\tIntersect_Count\tFile_A_Uni\tFile_B_Uni\tRho\tP\tRho_int\tP_int\n')

    unique_matrix = {}
    for file_a in file_list:

        # Restrict to max_peak_count lines
        max_file_a = add_name(file_a, pre='', suf=f'.{max_peak_count}', ext='')
        if not path.exists(max_file_a):
            !head -$max_peak_count $file_a > $max_file_a
        
        # Sort by chr, pos
        s_max_file_a = add_name(max_file_a, pre='', suf=f'.gsort', ext='')
        if not path.exists(s_max_file_a):
            !cat $max_file_a | sort -k1,1 -k2,2n > $s_max_file_a

        for file_b in file_list:
            
            # Add files to unique matrix  
            if file_a not in unique_matrix:
                unique_matrix[file_a] = {file_b: 0}
            elif file_b not in unique_matrix[file_a]:
                unique_matrix[file_a][file_b] = 0
            
            if file_b not in unique_matrix:
                unique_matrix[file_b] = {file_a: 0}
            elif file_a not in unique_matrix[file_b]:
                unique_matrix[file_b][file_a] = 0
                
            if unique_matrix[file_b][file_a] == 0:
            
                # Restrict to max_peak_count lines
                max_file_b = add_name(file_b, pre='', suf=f'.{max_peak_count}', ext='')
                if not path.exists(max_file_b):
                    !head -$max_peak_count $file_b > $max_file_b

                # Sort by chr, pos
                s_max_file_b = add_name(max_file_b, pre='', suf=f'.gsort', ext='')
                if not path.exists(s_max_file_b):
                    !cat $max_file_b | sort -k1,1 -k2,2n > $s_max_file_b

                # Only assess different files        
                if file_a == file_b:
                    pass
                else:

                    file_a_simp = file_a.split('/')[-1]
                    file_b_simp = file_b.split('/')[-1]
                    print(f'A: {file_a_simp}')
                    print(f'B: {file_b_simp}')

                    # Obtain peak ranks
                    a_peak_to_rank = {}
                    fa = open(max_file_a, 'r')
                    for a_rank, line in enumerate(fa):
                        a_peak_name = line.split('\t')[3]
                        a_peak_to_rank[a_peak_name] = a_rank
                    fa.close()

                    b_peak_to_rank = {}
                    fb = open(max_file_b, 'r')
                    for b_rank, line in enumerate(fb):
                        b_peak_name = line.split('\t')[3]
                        b_peak_to_rank[b_peak_name] = b_rank
                    fb.close()

                    # Perform AB intersections
                    file_ab_int = 'temp_AB_int.txt'
                    !bedtools intersect -b $s_max_file_a -a $s_max_file_b > $file_ab_int
                    file_ba_int = 'temp_BA_int.txt'
                    !bedtools intersect -b $s_max_file_b -a $s_max_file_a > $file_ba_int

                    # Open temp file
                    cor_file = 'temp_cor.txt'
                    corf = open(cor_file, 'w')

                    # Record intersection coords
                    # Coords are unique and same between AB and BA
                    # Each coord has one peak name
                    # A peak name can appear twice

                    coord_to_b_name = {}
                    ab_int = open(file_ab_int, 'r')
                    for line_count, line in enumerate(ab_int):
                        tab = line.split('\t')
                        peak_coord = f'{tab[0]}_{tab[1]}_{tab[2]}'
                        b_peak = tab[3]
                        coord_to_b_name[peak_coord] = b_peak
                    ab_int.close()

                    a_peak_set = set()
                    b_peak_set = set()
                    a_rank_list = []
                    b_rank_list = []
                    int_count = file_a_uni = file_b_uni = 0
                    ba_int = open(file_ba_int, 'r')
                    for line_count, line in enumerate(ba_int):
                        tab = line.split('\t')
                        peak_coord = f'{tab[0]}_{tab[1]}_{tab[2]}'
                        a_peak = tab[3]
                        b_peak = coord_to_b_name[peak_coord]

                        # Get the ranks
                        a_peak_set.add(a_peak)
                        b_peak_set.add(b_peak)
                        a_rank = a_peak_to_rank[a_peak]
                        b_rank = b_peak_to_rank[b_peak]
                        a_rank_list.append(a_rank)
                        b_rank_list.append(b_rank)
                        int_count += 1
                        corf.write(f'{a_peak}\t{a_rank}\t{b_rank}\t{b_peak}\n')

                    ba_int.close()
                    rho_int, p_int = spearmanr(a_rank_list, b_rank_list)
                    
                    # Add the peaks outside the intersection
                    for a_peak in a_peak_to_rank:
                        if a_peak not in a_peak_set:
                            file_a_uni += 1
                            a_rank = a_peak_to_rank[a_peak]
                            a_rank_list.append(a_rank)
                            b_rank_list.append(max_peak_count)
                            corf.write(f'{a_peak}\t{a_rank}\t{max_peak_count}\t.\n')

                    for b_peak in b_peak_to_rank:
                        if b_peak not in b_peak_set:
                            file_b_uni += 1
                            b_rank = b_peak_to_rank[b_peak]
                            a_rank_list.append(max_peak_count)
                            b_rank_list.append(b_rank)
                            corf.write(f'.\t{max_peak_count}\t{b_rank}\t{b_peak}\n')

    #                 corf.close()
                    rho, p = spearmanr(a_rank_list, b_rank_list)
                    print(rho, p)
                    out.write(f'{file_a}\t{file_b}\t{file_a_simp}\t{file_b_simp}\t{int_count}\t{file_a_uni}\t{file_b_uni}\t{rho}\t{p}\t{rho_int}\t{p_int}\n')
                    unique_matrix[file_a][file_b] = 1
                    unique_matrix[file_b][file_a] = 1

    out.close()  


def add_name(path_and_file_name, pre, suf, ext):
    """
    Rename and input file, respecting the path and extension

    Parameters
    ----------
    path_and_file_name : str
        path and file name, e.g., /path/to/file_name.txt
    pre : str
        prefix to be added before the file name
    suf : str
        suffix to be added after the file name
    ext : str
        extension to be added at the end of the file name

    Returns
    -------
    the_output : str
        modified path and file name, e.g., /path/to/prefix_file_name_suffix.extension

    """
    bits = path_and_file_name.split('/')
    path = '/'.join(bits[:-1])
    if len(path) > 0:
        path = f'{path}/'

    file_name_no_path = bits[-1]
    file_bits = file_name_no_path.split('.')

    recognized_ext = ['csv', 'txt', 'pdf', 'png', 'bed', 'fa']

    # Add a period to the extension as needed
    extension_mod = ''
    if ext is not None and ext != '':
        if ext.startswith('.'):
            extension_mod = ext
        else:
            extension_mod = f'.{ext}'

    the_output = ''
    if len(file_bits) == 1:
        # No periods in the name, therefore no extension
        if extension_mod is not None and extension_mod != '':
            the_output = f'{path}{pre}{file_name_no_path}{suf}{extension_mod}'
        else:
            the_output = f'{path}{pre}{file_name_no_path}{suf}'

    else:
        # Periods present in file name
        orig_ext = file_bits[-1]

        if orig_ext not in recognized_ext:
            # String after the period does not appear to be an extension
            if extension_mod is not None and extension_mod != '':
                the_output = f'{path}{pre}{file_name_no_path}{suf}{extension_mod}'
            else:
                the_output = f'{path}{pre}{file_name_no_path}{suf}'
        else:
            file_name_simp = '.'.join(file_bits[:-1])
            if extension_mod is not None and extension_mod != '':
                # New extension replaces the old one
                the_output = f'{path}{pre}{file_name_simp}{suf}{extension_mod}'
            else:
                # Keep the old extension
                the_output = f'{path}{pre}{file_name_simp}{suf}.{orig_ext}'

    if the_output != path_and_file_name:
        print(f'add_name output: {the_output}')
        return the_output
    else:
        print(f'ERR: input name and output name are the same.\n'
              f'Input: {path_and_file_name}\n'
              f'Output: {the_output}')
        return None      
        