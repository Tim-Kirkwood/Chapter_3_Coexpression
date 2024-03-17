# Chapter_3_Coexpression

## Data preparation
Count table - Go to Sulheim S et al. (2020).  Download Data S4, extract average values for all triplicate M145 datapoints for use as input to transcriptominer_thesis.py - note, input file is provided here as average_data.csv.  Place in same directory as this script.
Genome reference - Download NC_003888 from NCBI in genbank format, place in same directory as this script.
BGCs - download genbank files of BGCs of interest and place into bgc_files directory.

## Run program
Enter reference genome filename for map_file_name variable (line 31) and count_file_name (line 37).
Run transcriptominer_thesis.py (parameters between lines 16 and 50 can be changed according to preferences).  


Sulheim S, Kumelj T, van Dissel D, Salehzadeh-Yazdi A, Du C, van Wezel GP, Nieselt K, Almaas E, Wentzel A, Kerkhoven EJ. Enzyme-Constrained Models and Omics Analysis of Streptomyces coelicolor Reveal Metabolic Changes that Enhance Heterologous Production. iScience. 2020 Sep 3;23(9):101525. doi: 10.1016/j.isci.2020.101525. PMID: 32942174; PMCID: PMC7501462.
