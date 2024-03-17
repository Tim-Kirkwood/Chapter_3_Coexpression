# Chapter_3_Coexpression

## Data preparation
Count table - Go to Sulheim S et al. (2020).  Download Data S4, extract average values for all triplicate M145 datapoints for use as input to transcriptominer_thesis.py - note, input file is provided here as average_data.csv.  Place in same directory as this script.
Genome reference - Download NC_003888 from NCBI in genbank format, place in same directory as this script.
BGCs - download genbank files of BGCs of interest and place into bgc_files directory.

## Run program
Enter reference genome filename for map_file_name variable (line 31) and count_file_name (line 37). 
Run transcriptominer_thesis.py (extension settings and other parameters between lines 16 and 50 can be changed according to preferences).  

## Results
For each BGC in the bgc_files directory, 4 interactive html files will be written to the same directory as this script.  Filename formats shown below, bold substrings are dynamic and change depending on the query file/run parameters:
- **BGCfilename**_unfiltered_cor: Correlation matrix of neighborhood spanning from BGC start - ```upstream extension``` to BGC end - ```downstream extension``` (default extenions are 20 genes each way).  
- **BGCfilename**_mibig_cor: Correlation matrix of neighborhood spanning from BGC start to BGC end
- **BGCfilename**_line: Line plot of counts for genes in BGC ignoring extensions.  Either target (i) all genes, (ii) all genes of a certain functional annotation - see ```kinds``` in line 46, or (iii) user-specified genes - see ```user_line_cds``` in line 44.
- **BGCfilename**\_cor_**mask**: Correlation matrix of neighborhood spanning from BGC start - ```upstream extension``` to BGC end - ```downstream extension``` (default extenions are 20 genes each way).  Correlations below mask threshold are set to 0.  Mask threshold is **average pairwise correlation of bgc - noise** (mask is given in filename).  

Note - **BGCfilename** is the filename of a single bgc file from the bgc_files directory that is undergoing analysis.


Sulheim S, Kumelj T, van Dissel D, Salehzadeh-Yazdi A, Du C, van Wezel GP, Nieselt K, Almaas E, Wentzel A, Kerkhoven EJ. Enzyme-Constrained Models and Omics Analysis of Streptomyces coelicolor Reveal Metabolic Changes that Enhance Heterologous Production. iScience. 2020 Sep 3;23(9):101525. doi: 10.1016/j.isci.2020.101525. PMID: 32942174; PMCID: PMC7501462.
