# five_tr_chip
Code for analysis of ChIP-seq data from five transcription regulators

**Code for assessing whether peaks are in the proximal or distal region**
Example usage:

`python promoter_enhancer_split.py arid1b_peaks_hg38.txt human ARID1B`

**Code for assessing overlap between replicates**
Example usage:

`python compare_replicates.py -a arid1b_peaks_hg38.txt -b atac_peaks_hg38.txt -mo 0.01`

**Code for assessing and visualizing the overlap between ChIP peak sets**
Example usage:

`python analysis_code_NeuroDevEpi.py -m 0.4 -t hg38_GW23 -g ARID1B,BCL11A,FOXP1,TBR1,TCF7L2 -p arid1b_peaks_hg38.txt,bcl11a_peaks_hg38.txt,foxp1_peaks_hg38.txt,tbr1_peaks_hg38.txt,tcf7l2_peaks_hg38.txt -s 3.9`

**Function for the assessment of overlap between multiple files of ChIP-seq data (deployed in Jupyter Lab)**
Example usage:

`spearman_peaks(prom_liver_brain_list, 10000, 'prom_liver_brain_cor_10000.txt')`

Where:
'prom_liver_brain_list' is a python list of files of ChIP-seq peaks sorted by log10P 
'10000' is the number of peaks to use across the peak files
'prom_liver_brain_cor_10000.txt' is the output file name
