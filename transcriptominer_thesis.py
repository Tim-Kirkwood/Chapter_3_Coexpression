# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 11:08:32 2023

@author: u03132tk
"""
import time
import pandas as pd
import numpy as np
import os 
import copy
import plotly.express as px
from Bio import SeqIO

#what feature in the genome is used to identify proteins in the count table?    
flag = 'locus_tag'

#note this data was the averaged data from the publication Sulheim S, Kumelj T, van Dissel D, Salehzadeh-Yazdi A, Du C, van Wezel GP, Nieselt K, Almaas E, Wentzel A, Kerkhoven EJ. Enzyme-Constrained Models and Omics Analysis of Streptomyces coelicolor Reveal Metabolic Changes that Enhance Heterologous Production. iScience. 2020 Sep 3;23(9):101525. doi: 10.1016/j.isci.2020.101525. PMID: 32942174; PMCID: PMC7501462.
count_file = 'C:/Users/u03132tk/.spyder-py3/transcriptome_borders/average_data.csv'
extension_down = 20
extension_up = 20
noise = 0.1
user_line_cds = []
kinds = []# ['biosynthetic']
#genome accession NC_003888 - used in publication above as mapping reference for rna seq
map_file = rf"{os.getcwd()}\sequence (15).gb"
bgc_dir = rf"{os.getcwd()}\bgc_files"
x_axis = [21, 29, 33, 37, 41, 45, 49, 53, 57]

print ('processing data')
process_start = time.time()
map_records = list(SeqIO.parse(open(map_file,"r"), "genbank") )
bgc_files = os.listdir(bgc_dir)
for file in bgc_files:
    bgc_name = file[0:file.rindex('.')]
    print (bgc_name)
    bgc_file = f'{bgc_dir}/{file}'
    bgc_record = list(SeqIO.parse(open(bgc_file,"r"), "genbank") )
    assert len (bgc_record) == 1
    bgc_cds = [f.qualifiers['translation'] for f in bgc_record[0].features if f.type == 'CDS']
    line_cds = []
    print ('making line cds')
    for f in bgc_record[0].features:
        if f.type == 'CDS':
            if 'gene_kind' in f.qualifiers.keys():
                if kinds != []:
                    
                    if any([i in kinds for i in f.qualifiers['gene_kind']]):
                        line_cds += [f.qualifiers['translation']]
                else:
                    
                    line_cds += [f.qualifiers['translation']]
    
    #match = []
    for record in map_records:
        #assumes that the order of features matches genomic order
        map_cds = [(f.qualifiers[flag], f.qualifiers['translation']) for f in record.features if f.type == 'CDS']
        match = [(i, v[0]) for i, v in enumerate(map_cds) if v[1] in bgc_cds]
        groups = []
        
        group = []
        #can have multiple copies of a gene
        #so split up all of the ones witth match proteins that arent adjacent and only keep ones with same size as bgc
        for index, entry in enumerate(match):
            if index != len(match) - 1:
                #print (index, abs(match[index + 1][0] - entry[0]))
                if abs(match[index + 1][0] - entry[0]) > 1:
                    #print ('true')
                    #it is close enough to previous entry so add to current group
                    group += [entry]
                    #this group is done so finish
                    groups += [copy.deepcopy(group)]
                    #re-initiate
                    group = []
                    
                else:
                    #print ('false')
                    group += [entry]
                    #print (group)
            else:
                if group == []:
                    #the last group got seperated so add alone
                    groups += [[entry]]
                else:
                    #it is close enough to previous entry so add to current group
                    group += [entry]
                    #this group is done so finish
                    groups += [copy.deepcopy(group)]
        potential_groups = []
        for group in groups:
            if len(group) == len(bgc_cds):
                potential_groups += [group]
        assert len(potential_groups) == 1
        #raise ValueError
        match = potential_groups[0]    
        locii = [i[0] for i in match]
        
        
        span = max(locii) - min(locii) + 1 #0-indexed
        assert span == len(bgc_cds), f"span {span} bgc {len(bgc_cds)}"#they are adjacent
        mibig_window = [i[0][0] for i in map_cds[min(locii) : max(locii)+1]]
        transcriptome_window = [i[0] for i in map_cds[min(locii) - extension_up : max(locii) + extension_down]]
        assert all([len(i) == 1 for i in transcriptome_window])
        transcriptome_window = [i[0] for i in transcriptome_window]
        line_locii = [i[0] for i in match if map_cds[i[0]][1] in line_cds]
        line_cds = [map_cds[l][0][0] for l in line_locii]
        #raise ValueError
        break
    assert match != []
    
    
    count_table = pd.read_csv(count_file, 
                              sep = ',',
                              header = 0,
                              index_col = 0).transpose()
    assert len(set(count_table.columns)) == len(count_table.columns) 
    bgc_table = count_table[transcriptome_window]
    #plot lines so you can identify trends
    
    pearson_start = time.time()
    #might need to normalise by gene length!
    cor_table = bgc_table.corr().fillna(0)#correlationg 2 0 arrays give nan
    fig = px.imshow(cor_table)
    fig.write_html(f'{bgc_name}_unfiltered_cor.html')
    mibig_table = count_table[mibig_window]
    if user_line_cds != [] or line_cds != []:
        print ('writing line graph')
        line_df = []
        if user_line_cds != []:
            reference = user_line_cds
        else:
            reference = line_cds
        for gene in mibig_table.columns:
            if gene in reference:
                y = mibig_table[gene].values.tolist()
                assert len(y) == len(x_axis)
                line_df += [[gene, x, y] for x, y in zip(x_axis, y)]
        line_df = pd.DataFrame(line_df, columns = ['gene', 'time (h)', 'normalised count'])
        fig = px.line(line_df, x="time (h)", y="normalised count", color='gene')
        fig.write_html(f'{bgc_name}_line.html')
    
    
    mibig_cors = mibig_table.corr().fillna(0)
    fig = px.imshow(mibig_cors)
    fig.write_html(f'{bgc_name}_mibig_cor.html')
    mibig_mask = mibig_cors.abs().mean().mean() - noise
    print ('mask:  ', mibig_mask)
    cor_table = cor_table.mask((cor_table < mibig_mask) & (cor_table > -mibig_mask)).fillna(0)
    signs = np.sign(cor_table)
    fig = px.imshow(cor_table)
    fig.write_html(f'{bgc_name}_cor_{mibig_mask}.html')