#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:58:54 2021

@author: mheinzinger
"""

from pathlib import Path

# for tensor operations, data loading etc
import numpy as np
import h5py
from pandas import DataFrame
from scipy.spatial.distance import pdist, squareform

# for dimensionality reduction
import umap

# for visualzation
import matplotlib.pyplot as plt
import seaborn as sns # only necessary if you want to play around with plots
import plotly
import csv


def save_plotly_figure_to_html(figure, path):
    """
    Store plotly figure as interactive HTML file
    :param figure: A Plotly Figure
    :param path: A string representing the path and/or filename where the HTML figure should be stored
         (e.g.: /path/to/figure.html).
    """
    plotly.offline.plot(figure, filename=path)
    return None

# https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
def render_3D_scatter_plotly(embeddings_dataframe):
    """
    Return a Plotly Figure (3D scatter plot) based on a DataFrame containing three components.
    :param embeddings_dataframe: the DataFrame *must* contain three numerical columns called `component_0`,
            `component_1` and `component_2`. The DataFrame index will be used to identify the points in the
            scatter plot. Optionally, the DataFrame may contain a column called `label` which will be used
            to color the points in the scatter plot.
    :return: A 3D scatter plot
    """
    import plotly.express as px

    fig = px.scatter_3d(embeddings_dataframe,
                            x='component_0',
                            y='component_1',
                            z='component_2',
                            color='label',
                            symbol='label',
                            hover_name=embeddings_dataframe.index,
                            hover_data=["hover_data"],
                            # This tries to avoid repetition of symbols:
                            symbol_sequence = ['circle', 'square', 'diamond', 
                                               'cross', 'x', 'circle-open', 
                                               'square-open', 'diamond-open', 
                                               ],
                            )

    fig.update_layout(
        # Remove axes ticks and labels as they are usually not informative
        scene=dict(
            xaxis=dict(
                showticklabels=False,
                showspikes=False,
                title=""
            ),
            yaxis=dict(
                showticklabels=False,
                showspikes=False,
                title=""
            ),
            zaxis=dict(
                showticklabels=False,
                showspikes=False,
                title=""
            )
        ),
    )

    return fig

def pairwise_distances(data, metric="euclidean"):
    # usually euclidean or cosine distance worked best
    return pdist(data,metric=metric)

def read_csv(label_csv_p, representatives, id_idx=0):
    char10_to_labels = dict()
    with open(label_csv_p,'r') as in_f: # open CSV
        next(in_f) # skip header: 10char;FASTA;Name;Major group;Species;Major taxon
        reader = csv.reader(in_f,delimiter=",") # read in CSV
        for line in reader:
            identifier = line[id_idx]
            label = line[4] # major_group=4; major_taxon=6
            char10_to_labels[identifier] = (label," ++ ".join(line[2:]))
    return char10_to_labels


def read_embeddings( emb_p):
    # load pre-computed embeddings in .h5 file format
    # returns dictionary with fasta headers as keys and a single vector (embeddings) per protein as values
    # values have 1024-dimensions for ProtT5 and 128-dimensions for ProtTucker
    print("Loading pre-computed embeddings from: {}".format(emb_p))
    h5_f = h5py.File(emb_p,'r')
    dataset = { pdb_id : np.array(embd) for pdb_id, embd in h5_f.items() }
    print("Example: {}".format(next(iter(dataset.keys()))))
    print("Number of embeddings: {}".format(len(dataset)))
    return dataset

def read_fasta_ids(rep_seqs):
    with open(rep_seqs,'r') as in_f:
        # replace weird characters that mess up H5 format (required for mapping to embeddings)
        sequence_headers = { line.strip().replace(">","").replace("/","_").replace(".","_") for line in in_f if line.startswith(">") }
    print("Read {} sequences from FASTA.".format(len(sequence_headers)))
    return sequence_headers

def read_id_mapping(id_mapping_p):
    with open(id_mapping_p,"r") as in_f:
        fastaID_2_char10 = { line.strip().split()[1] : line.strip().split()[0] for line in in_f }
    return fastaID_2_char10

def main():
    # root directory that holds, proteins.fasta, embeddings.h5, labels.csv and some output_file.html
    root = Path.cwd() / "mysite"
    rep_seqs = root / "3and6_10char.fasta"
    emb_p = root / "3and6_10char.h5" 
    label_csv_p = root / "3and6_w10Char.csv"

    fig_3D_p = str(root / "3and6_3D_majorTaxon.html")
    fig_2D_p = str(root / "3and6_2D_majorTaxon.pdf")

    # read in proteins.fasta to get IDs of all proteins
    representatives = read_fasta_ids(rep_seqs)
    # read in labels.csv to retrieve annotations for all proteins in representatives
    grouping = read_csv(label_csv_p, representatives, id_idx=0)

    raw_embeddings = read_embeddings( emb_p) # reads in embeddings from H5PY format
    embeddings = { identifier : embd if identifier in grouping else print(identifier) for identifier, embd in raw_embeddings.items() }
    
    keys, data = zip(*embeddings.items())
    data = np.vstack(data) # matrix of values (protein-embeddings); n_proteins x embedding_dim

    #data = np.random.randn(*data.shape) # uncomment this line if you want to replace embeddings with random vectors
    
    # data should be n_proteins x 1024 (ProtT5) OR n_proteins x 128 (ProtTucker)
    print("Shape of raw embeddings (num_proteins x embedding dimension): {}".format(data.shape))
    
    # get pairwise distances; should be n_proteins x n_proteins
    pdist = squareform( pairwise_distances(data) )
    print("Shape of pairwise distance matrix (num_proteins x num_proteins): {}".format(pdist.shape)) 
    
    # visualize high-dimensional embeddings with dimensionality reduction (here: umap)
    # Tutorial: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
    # Parameters: https://umap-learn.readthedocs.io/en/latest/parameters.html
    fit = umap.UMAP(n_neighbors=25, min_dist=0.5, random_state=42, n_components=3) # initialize umap; use random_state=42 for reproducability
    u = fit.fit_transform(data) # fit umap to our embeddings
    
    
    raw_data = dict()
    labels = list()
    hover_data = list()
    for idx, key in enumerate(keys):
        group, more_info = grouping[key]
        if group not in raw_data:
            raw_data[group] = list()
        raw_data[group].append(u[idx,:2]) # Only necessary for 2D plot
        hover_data.append(more_info)
        labels.append(group)
        
        
    embeddings_dataframe = DataFrame(u, columns=["component_0", "component_1", "component_2"])
    embeddings_dataframe['label'] = labels
    embeddings_dataframe['hover_data'] = hover_data
    #embeddings_dataframe.index = data.index
    fig = render_3D_scatter_plotly(embeddings_dataframe=embeddings_dataframe)
    
    fig.update_traces(marker=dict(size=6,
                                  line=dict(width=1,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))

    
    save_plotly_figure_to_html(fig, fig_3D_p)
    
    fit = umap.UMAP(n_neighbors=25, min_dist=0.5, random_state=42, n_components=2) # initialize umap; use random_state=42 for reproducability
    u = fit.fit_transform(data) # fit umap to our embeddings
    
    raw_data = dict()
    for idx, key in enumerate(keys):
        group, _ = grouping[key]
        if group not in raw_data:
            raw_data[group] = list()
        raw_data[group].append(u[idx,:2]) # Only necessary for 2D plot
    
    
    for group, umap_data in raw_data.items():
        raw_data[group] = np.vstack(umap_data)
        ax = sns.scatterplot(x=raw_data[group][:,0], y=raw_data[group][:,1], label=group)
    
    
    plt.legend(prop={'size': 10})
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.figure.savefig(fig_2D_p, bbox_inches='tight')
    

if __name__ == '__main__':
    main()