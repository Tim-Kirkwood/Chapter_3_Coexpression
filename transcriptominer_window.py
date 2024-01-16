# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 11:08:32 2023

@author: u03132tk
"""
import time
import pandas as pd
import numpy as np
import networkx as nx
import itertools
import plotly.graph_objects as go
from random import randint
import sys
import statistics
import os 
import copy
import plotly.express as px
import re
from Bio import SeqIO
import copy
import numpy as np

def weighted_j_sim(array1, array2):
    return np.minimum(array1, array2).sum()/np.maximum(array1, array2).sum()
def make_wjs_df(df):
    wjs_matrix = []
    for index, (name, values) in enumerate(df.iterrows()):
        row = []
        for other_index, (other_name, other_values) in enumerate(df.iterrows()):
            #if other_index>index: #dont check self or something already compared
            weighted_js = weighted_j_sim(values, 
                                             other_values)
            row+= [weighted_js]
        wjs_matrix += [row]
    return pd.DataFrame(wjs_matrix, columns = df.columns, index = df.index)
    
def build_edge_traces(G : nx.classes.graph.Graph, df_all, weight_label = '# of shared modules'): 
    all_weights = [weight for weight in set(df_all['weight']) if weight != 0]
    if all_weights == []:
        return []
    edge_traces = []
    # try:
    #     dcolorscale=[[weight/max(all_weights), f"rgb({255*(weight/max(all_weights))},0,0)"] 
    #              for weight in range(0, max(all_weights) + 1)]
    # except TypeError as e:
    #     print (locals())
    #     raise e
    for weight in all_weights:
        df = df_all[df_all['weight'] == weight]
        weight = set(df['weight'])
        #assert len(weight) == 1
        weight = list(set(weight))[0]
        #bool_index = unique_edges['weight'] == weight
        #weight_df = unique_edges[bool_index]
        edge_x = []
        edge_y = []
        middle_x = []
        middle_y = []
        middle_text = []
        for _, edge in df.iterrows():
            x0, y0 = G.nodes[edge[0]]['pos']
            x1, y1 = G.nodes[edge[1]]['pos']
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
            middle_x += [(x0+x1)/2]
            middle_y += [(y0+y1)/2]
            middle_text += [f'{weight_label} = {weight}']
        #TODO edge color bar
        edge_traces.append(go.Scatter(
                                        x=middle_x,
                                        y=middle_y,
                                        text=middle_text,
                                        mode='markers',
                                        hoverinfo='text',
                                        marker=go.Marker(
                                            opacity=0
                                            )
                                        ))
        edge_traces.append(go.Scatter(x=edge_x, y=edge_y,
                                      #scaling val is number of prots for which this multiplication works well
                                      line=dict(width = 0.1+(3*weight/max(all_weights)),#*scaling_val/len(formatted_unique_edges), 
                                                color = f"rgba({round(255*(weight/max(all_weights)))},0,0,1)"),#{0.6 + (0.4 * (weight/max(all_weights)))}
                                      #marker=dict(#showscale=True,
                                                  
                                                  #colorscale=dcolorscale,
                                                  
                                       #           size=10000#,
                                                  # colorbar=dict(thickness=15,
                                                  #               title='Edge weights',
                                                  #               xanchor='left',
                                                  #               titleside='right'
                                                  #               )
                                         #         ),
                                      
                                      mode='lines'))#https://stackoverflow.com/a/61349739/11357695#
    return edge_traces

def build_node_trace(G : nx.classes.graph.Graph):
    #TODO add option to colour by protein annotation
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)
    return go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))
def node_trace_format_data (G : nx.classes.graph.Graph):
    annotations = []
    for node in G.nodes():
        annotations+=[node]
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))#will be bad if one prot has manyyyy more connections than pothers
        node_text.append(f'{annotations[node]} - # of connections: {str(len(adjacencies[1]))}')
    return {'node_adjacencies' : node_adjacencies, 
            'node_text' : node_text}

def plot_graph(edge_table, filepath, write_to_file = True):
    G = nx.from_pandas_edgelist(edge_table, 
                                edge_attr=True) 
    G = G.to_undirected()
    pos = nx.drawing.layout.spring_layout(G)
    nx.set_node_attributes(G, pos, 'pos')
    if not write_to_file or filepath == None:   
        return G
    edge_traces = build_edge_traces(G, 
                                       edge_table)
    node_trace = build_node_trace(G)
    
    
    
    formatting_data = node_trace_format_data(G)
    node_trace.marker.color = formatting_data['node_adjacencies']
    node_trace.text = formatting_data['node_text']
    
    fig = go.Figure(data=edge_traces + [node_trace],
                 layout=go.Layout(
                    title='<br>Network graph made with Python',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                    #    text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.write_html(filepath)

def axes_to_gap_matrix(axes):
    matrix = []
    for index, label in enumerate(axes):
        end_val = len(axes) - index
        start_val = -index
        matrix += [list(range(start_val, end_val))]
    return matrix






def check_tuple_list_overlap(tuple_indexes: list):
    for check_index, tup_check in enumerate(tuple_indexes):
        #print (f'\n\ntup_check: {tup_check}\n')
        start_check, end_check = tup_check
        for compare_index, tup_compare in enumerate(tuple_indexes):
            if compare_index != check_index:
                start_compare, end_compare = tup_compare
                if start_compare <= start_check <= end_compare:
                    return True
                if start_compare <= end_check <= end_compare:
                    return True
                if start_compare >= start_check and end_compare <= end_check:
                    return True
                if start_compare <= start_check and end_compare >= end_check:
                    return True
    return False



def location_list(string_location, need_strand=True):
    '''
    split a json/gbk location id into list of index ints - e.g. location id:  join{[22623:>22761](-), [22289:22430](-), [22155:22299](-)}
    '''
    #TODO - check if the strand affects the locii - should not do, biopython at least keeps with + strand index
    obj = re.split('[^0-9]', string_location)
    obj_filtered = [item for item in obj if item != '']
    obj_int = list(map(int, obj_filtered))
    feature_type = 'good'
    if 'join' in string_location:
        # join indexes will have intermediate indexes as well, which are ignored in terms of outputting total span
        obj_int = [min(obj_int), max(obj_int)]
        feature_type = 'join'
    else:
        if len(obj_filtered) > 2:
            raise ValueError('feature has >2 parts but does not contain a "join" flag - unable to confirm whether this feature should be treated as a single or multi-part span')
    return [feature_type, obj_int]

# =============================================================================
# ---BUGS---
# =============================================================================
#TODO - some BGCs have matrixes where there are no vales > min core for last coupl eof columns/rows - eg bgc EV45_RS29715.  Also this bgc should include the protein after the terminal 0 but doesnt?
#seeems like an indexing error when you went table[min(proteins):max(proteins)] - should be max(proteins)+1

flag = 'locus_tag'
extension_down = 20
extension_up = 20
noise = 0.1
x_axis = [21, 29, 33, 37, 41, 45, 49, 53, 57]
user_line_cds = []# ["SCO2776", "SCO2777", "SCO2778", "SCO2779", "SCO2780", "SCO2781", "SCO2782", "SCO2783", "SCO2784", "SCO2785", "SCO2786"]#"SCO7461",
# "SCO7462",
# "SCO7463",
# "SCO7464",
# "SCO7465",
# "SCO7466",
# "SCO7467",
# "SCO7468",
# "SCO7469",
# "SCO7470",
# "SCO7471",
# "SCO7472",
# "SCO7473",
# "SCO7474",
# "SCO7475"]
kinds = []# ['biosynthetic']
print ('processing data')
process_start = time.time()
map_file = "C:/Users/u03132tk/.spyder-py3/transcriptome_borders/sequence (15).gb"
map_records = list(SeqIO.parse(open(map_file,"r"), "genbank") )

bgc_dir = "C:/Users/u03132tk/.spyder-py3/transcriptome_borders/bgc_files"
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
                    print ('adding all')
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
    
    count_file = 'C:/Users/u03132tk/.spyder-py3/transcriptome_borders/average_data.csv'
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
    
    cor_wjs_df = make_wjs_df(cor_table.abs())
    fig = px.imshow(cor_wjs_df * signs)
    fig.write_html(f'{bgc_name}_cor_wjs.html')
    
    #cor_table = cor_table.mask(cor_table < 0.7).fillna(0)
    r2_table = pd.DataFrame(cor_table.values * cor_table.values,
                            columns= cor_table.columns,
                            index = cor_table.index)
    fig = px.imshow(r2_table * signs)
    fig.write_html(f'{bgc_name}_r2.html')
    
    r2_wjs_df = make_wjs_df(r2_table)
    fig = px.imshow(r2_wjs_df * signs)
    fig.write_html(f'{bgc_name}_r2_wjs.html')
    #print (cor_table.head())

    