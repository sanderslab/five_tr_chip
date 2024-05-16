# five_tr_chip
Code for analysis of ChIP-seq data from five transcription regulators



Function for the assessment of overlap between multiple files of ChIP-seq data (deployed in Jupyter Lab)
Example usage:
<spearman_peaks(prom_liver_brain_list, 10000, 'prom_liver_brain_cor_10000.txt')>

Where:
'prom_liver_brain_list' is a python list of files of ChIP-seq peaks sorted by log10P 
'10000' is the number of peaks to use across the peak files
'prom_liver_brain_cor_10000.txt' is the output file name
