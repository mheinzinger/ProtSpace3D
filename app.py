#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 23:00:13 2022

@author: mheinzinger
"""

from pathlib import Path
from collections import defaultdict
import base64
import re

# for tensor operations, data loading etc
import numpy as np
import h5py
from pandas import DataFrame

# for dimensionality reduction
import umap
import csv

# for visualzation
import dash
from jupyter_dash import JupyterDash
from dash import dcc, html, Input, Output, State, no_update
import dash_bio as dashbio

import plotly.graph_objects as go
import dash_bio.utils.ngl_parser as ngl_parser
from dash.exceptions import PreventUpdate

symbols=['circle', 'square', 'diamond', 'cross', 'x', 
         'circle-open', 'square-open', 'diamond-open', 
        ]

class StructureContainer(object):
    def __init__(self, pdb_d):
        self.pdb_d = pdb_d

    def __call__(self):
        return self.public_seq_id

    def set_focus_point(self, curveNumber, pointNumber):
        self.curveNumber = curveNumber
        self.pointNumber = pointNumber
        return None

    def get_focus_point(self):
        return (self.curveNumber, self.pointNumber)

    def get_structur_dir(self):
        return self.pdb_d

    def set_structure_ids(self, seq_ids):
        #parser = PdbParser(pdb_p)
        #structure = parser.mol3d_data()
        if isinstance(seq_ids, list):
            self.public_seq_id = [seq_id.replace(
                ".", "_", 1) for seq_id in seq_ids]
        else:
            self.public_seq_id = [seq_ids.replace(".", "_", 1)]
        return None


# https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
def get_figure(df, filter_by, highlight=None):
    """
    Return a Plotly Figure (3D scatter plot) based on a DataFrame containing three components.
    :param embeddings_dataframe: the DataFrame *must* contain three numerical columns called `x`,
            `y` and `z`. The DataFrame index will be used to identify the points in the
            scatter plot. Optionally, the DataFrame may contain a column called `filter_by` which will be used
            to color the points in the scatter plot.
    :return: A 3D scatter plot
    """
    groups = set(df[filter_by])

    df["class_index"] = np.ones(len(df["SHORT_ID"]))*-100

    data = []
    for class_idx, group in enumerate(groups):
        df_group = df[df[filter_by] == group]
        trace = go.Scatter3d(x=df_group['x'],
                             y=df_group['y'],
                             z=df_group['z'],
                             mode='markers',
                             name=group,
                             # 10 colors are available; once those are used, pick different symbol
                             marker = dict(symbol = symbols[class_idx//10])
                             )
                             
        data.append(trace)
        df["class_index"][df[filter_by] == group]=class_idx
    if False:
        click_x, click_y, click_z = highlight
        trace = go.Scatter3d(x=[click_x], y=[click_y], z=[click_z],
                             mode='markers',
                             name="Focus",
                             )
        data.append(trace)

    fig = go.Figure(data=data)

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


def read_funfam_mapping(file_in):

    funfam_mapping = dict()

    with open(file_in) as read_in:
        next(read_in)
        for line in read_in:
            splitted_line = line.strip().split()
            raw_id = splitted_line[1]
            #uni_id = '{}/{}'.format(splitted_line[1], splitted_line[2])
            funfam_id = "FunFam: {} ## UniProt: {} ({})".format(
                splitted_line[0], splitted_line[1], splitted_line[2])
            funfam_mapping[raw_id] = funfam_id

    return funfam_mapping


def read_embeddings( emb_p):
    # load pre-computed embeddings in .h5 file format
    # returns dictionary with fasta headers as keys and a single vector (embeddings) per protein as values
    # values have 1024-dimensions for (T5_XL_U50_Pla2g2_alignment) and 128-dimensions for (PT5Tucker_Pla2g2_alignment)
    print("Loading pre-computed embeddings from: {}".format(emb_p))
    h5_f = h5py.File(emb_p,'r')
    dataset = { pdb_id : np.array(embd) for pdb_id, embd in h5_f.items() }
    print("Example: {}".format(next(iter(dataset.keys()))))
    print("Number of embeddings: {}".format(len(dataset)))
    return dataset

def read_csv(seq_ids,grouping_p):
    annotations = defaultdict(dict)
    add_info = dict()
    with open(grouping_p,'r') as read_in:
        # new_10_Char,old_10_char,Original fasta header,Name ,MajorGroup,Species,MajorTaxon
        next(read_in)
        reader = csv.reader(read_in)
        for splitted_line in reader:
            char10_id = splitted_line[0]
            if char10_id not in seq_ids:
                continue
            annotations[char10_id] = {
                'major_group': splitted_line[4],
                'species': splitted_line[5],
                'major_taxon': splitted_line[6],
            }
            add_info[char10_id] = "\n".join(splitted_line[2:])
    return annotations, add_info

def read_fasta_ids(rep_seqs):
    with open(rep_seqs,'r') as in_f:
        # replace special characters that could crash H5 loading
        sequence_headers = { line.strip().replace(">","").replace("/","_").replace(".","_") for line in in_f if line.startswith(">") }
    print("Read {} sequences from FASTA.".format(len(sequence_headers)))
    return sequence_headers

def read_data():
    root_dir = Path.cwd() / "mysite"
    rep_seqs = root_dir / "3and6_10char.fasta"
    emb_p = root_dir / "3and6_10char.h5"
    png_d = root_dir / "rank_1_10char_pngs"
    pdb_d = root_dir / "rank_1_10char"
    csv_p = root_dir / "3and6_w10Char.csv"

    seq_ids = read_fasta_ids(rep_seqs)
    structContainer = StructureContainer(pdb_d)
    grouping, additional_information = read_csv(seq_ids,csv_p)

    embeddings = read_embeddings( emb_p) # reads in embeddings from H5PY format
    embeddings = { identifier : embd if identifier in grouping else print(identifier)
                  for identifier, embd in embeddings.items()
                      }

    embeddings = read_embeddings(emb_p)  # reads in embeddings from H5PY format

    print("Number of embeddings before filtering for available grouping: {}".format(
        len(embeddings)))
    embeddings = {identifier: embd for identifier,
                  embd in embeddings.items() if identifier in grouping}
    print("Number of embeddings after filtering for available grouping: {}".format(
        len(embeddings)))

    keys, data = zip(*embeddings.items())
    # matrix of values (protein-embeddings); n_proteins x embedding_dim
    data = np.vstack(data)

    # data should be n_proteins x 1024 (PT5Tucker_Pla2g2_alignment) OR n_proteins x 128 (PT5Tucker_Pla2g2_alignment)
    print("Shape of raw embeddings (num_proteins x embedding dimension): {}".format(
        data.shape))

    # visualize high-dimensional embeddings with dimensionality reduction (here: umap)
    # Tutorial: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
    # Parameters: https://umap-learn.readthedocs.io/en/latest/parameters.html
    # initialize umap; use random_state=42 for reproducability
    # Initial attempt with n_neighbors=10, min_dist=0.3
    fit = umap.UMAP(n_neighbors=25, min_dist=0.5, 
                    random_state=42, n_components=3)
    u = fit.fit_transform(data)  # fit umap to our embeddings

    major_group = list()
    species = list()
    major_taxon = list()
    short_ids = list()
    hover_data=list()
    imgs = list()

    for idx, key in enumerate(keys):
        if idx == 0:
            print("Example identifier: {}".format(key))
        add_info = additional_information[key]
        hover_data.append(add_info)
        species.append(grouping[key]["species"])
        major_group.append(grouping[key]["major_group"])
        major_taxon.append(grouping[key]["major_taxon"])
        short_id = key.split("|")[1] if "|" in key else key.split()[0].replace(">", "")
        short_ids.append(short_id)
        png_p = png_d / (short_id+'.png')
        imgs.append("data:image/png;base64, " +
                    base64.b64encode(open(png_p, 'rb').read()).decode('ascii'))

    df = DataFrame(u, columns=["x", "y", "z"])

    df['hover_data'] = hover_data
    df['IMG_URL'] = imgs
    df['SHORT_ID'] = short_ids
    df['SPECIES'] = species
    df["MAJOR_GROUP"] = major_group
    df["MAJOR_TAXON"] = major_taxon
    fig = get_figure(df, filter_by="MAJOR_GROUP")
    return fig, df, structContainer


fig, df, structContainer = read_data()

# turn off native plotly.js hover effects - make sure to use
# hoverinfo="none" rather than "skip" which also halts events.
fig.update_traces(hoverinfo="none", hovertemplate=None)
fig.update_layout(clickmode='event+select', hovermode= 'closest')
fig.update_traces(marker=dict(size=6,
                              line=dict(width=1,
                                        color='DarkSlateGrey'),
                              ),
                  selector=dict(mode='markers'),
                  )

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = JupyterDash(__name__, external_stylesheets=external_stylesheets)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll',
        'overflowY': 'scroll',
        'width': '90vh',
        'height': '90vh'
    }
}

representation_options = [
    {"label": "backbone", "value": "backbone"},
    {"label": "ball+stick", "value": "ball+stick"},
    {"label": "cartoon", "value": "cartoon"},
    {"label": "hyperball", "value": "hyperball"},
    {"label": "licorice", "value": "licorice"},
    {"label": "axes+box", "value": "axes+box"},
    {"label": "helixorient", "value": "helixorient"},
    {"label": "rocket", "value":"rocket"},
    {"label": "rope", "value":"rope"},
    {"label": "surface", "value":"surface"},
    {"label": "tube", "value":"tube"},
    {"label": "distance", "value":"distance"},
    {"label": "contact", "value":"contact"},
    {"label": "helixorient", "value":"helixorient"},
]


dropdown_options = [{"label": seq_id, "value": seq_id}
                    for seq_id in df['SHORT_ID']]

app.layout = html.Div([
    html.Div([
        html.Div([
            dcc.Markdown("""
                    **3D scatter plot of embedding space.**
                """),
            dcc.Dropdown(
                ["SPECIES", "MAJOR_GROUP", "MAJOR_TAXON"],
                'MAJOR_GROUP',
                id='crossfilter-xaxis-column',
            ),
            dcc.Graph(
                id="graph",
                figure=fig,
                clear_on_unhover=True,
                style=styles['pre']
            ),

        ], className="six columns"),
        dcc.Tooltip(
            id="graph-tooltip"
        ),
        html.Div([
            html.Div([

                dcc.Markdown("""
                    **Click on points in the graph to show 3D structure here.**
                """),

                # dcc.Markdown("""
                #    **Choose your preferred representation:**
                # """),
                dcc.Dropdown(
                    id="nglstyle-dropdown",
                    options=representation_options,
                    multi=True,
                    value=["cartoon", "axes+box"]),

                # dcc.Markdown("""
                #    **Only relevenat if showing more than one protein. Allows you to avoid clashes:**
                # """),
                dcc.RadioItems(
                    id="nglstyle-radio",
                    options=[
                        {'label': 'sideByside', 'value': "True"},
                        {'label': 'Independent', 'value': "False"},
                    ],
                    value="False"
                ),
                # dcc.Markdown("""
                #    **Adjust the height and width of the window showing the 3D structure:**
                # """),
                dcc.Slider(
                    id='height-ngl-h',
                    min=300,
                    max=1000,
                    value=600,
                    step=100,
                    marks={300: '300px', 1000: '1000px'}
                ),
                dcc.Slider(
                    id='width-ngl-w',
                    min=300,
                    max=1000,
                    value=600,
                    step=100,
                    marks={300: '300px', 1000: '1000px'}
                ),
                # dcc.Markdown("""
                #    **Adjust quality of the 3D structure rendering:**
                # """),
                dcc.Dropdown(
                    id="ngl-stage-quality-dropdown",
                    value='auto',
                    options=[
                        {"label": s.capitalize(), "value": s}
                        for s in ["auto", "low", "medium", "high"]
                    ]
                ),
                # dcc.Markdown("""
                #    **Choose your protein(s) here:**
                # """),
                dcc.Dropdown(
                    id="default-ngl-molecule-dropdown",
                    options=dropdown_options,
                    placeholder="Select a molecule",
                    # value=dropdown_options[0]["value"],
                    multi=True
                ),

                # dcc.Markdown("""
                #             **Choose specific residues here (range has to start with ":");
                #             highlighted residues need to start with "@" and need to be separated by commas.**
                #    """),

                dcc.Input(
                    id="chain-atom-input",
                    placeholder="Eg. :629-819@700,750,800",
                    value=""
                ),
                dcc.Dropdown(
                    id="chain-atom-color",
                    placeholder="Select a color for highlighting specific residues.",
                    options=[{"label": s.capitalize(), "value": s} for s in ["black", "white", "red", "blue"]]),

                # , style=styles['pre']
                dashbio.NglMoleculeViewer(id="default-ngl-molecule"),

            ], className='six columns', style=styles['pre']),
        ], className="six columns"),
    ], className="row")
])


@app.callback(
    Output('graph', 'figure'),
    Input('crossfilter-xaxis-column', 'value'),
    Input("graph", "clickData"),  # event triggered by click on 3D graph
    # event triggered by dropdown
    Input("default-ngl-molecule-dropdown", "value"),
    State('graph', 'figure'),
)
def update_graph(xaxis_column_name, clickData, dropdownVal, fig):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    # If a click or dropdown action was triggered
    if ( ctx.triggered[0]["prop_id"] == "graph.clickData" or
        (ctx.triggered[0]["prop_id"] == "default-ngl-molecule-dropdown.value" and len(dropdownVal)>0 )):

        # if action was triggered by dropdown
        if ctx.triggered[0]["prop_id"] == "default-ngl-molecule-dropdown.value":
            seq_id=dropdownVal[0]
            entry = df[df["SHORT_ID"] == seq_id]
            curveNumber=int(entry["class_index"])
            print(curveNumber)
            class_points=df[df["class_index"]==curveNumber].reset_index()
            pointNumber=class_points.index[class_points["SHORT_ID"] == seq_id].tolist()[0]
            #pointNumber=df.index[(df["class_index"]==curveNumber & df["SHORT_ID"] == seq_id)]
            print(pointNumber)
            structContainer.set_focus_point(curveNumber,
                                            pointNumber)
            click_x, click_y, click_z = (float(entry["x"]),
                                         float(entry["y"]),
                                         float(entry["z"])
                                         )

        # if action was triggered by click --> get x,y,z of point
        elif ctx.triggered[0]["prop_id"] == "graph.clickData":
            structContainer.set_focus_point(clickData['points'][0]["curveNumber"],
                                            clickData['points'][0]["pointNumber"])
            click_x, click_y, click_z = (clickData['points'][0]['x'],
                                         clickData['points'][0]['y'],
                                         clickData['points'][0]['z']
                                         )

        # if there was already a highlighted/clicked point before: remove again
        if len(fig["data"])>len(np.unique(df["class_index"])):
            fig["data"] = [ fig["data"][i] for i in range(0, len(fig["data"])-1) ]
        # generate again figure from raw states
        fig = go.Figure(fig)
        # add new trace for highlighted/clicked point
        fig = fig.add_trace( go.Scatter3d(
                                    x=[click_x], y=[click_y], z=[click_z],
                                    mode='markers',
                                    name="Focus",
                                    marker=dict(
                                        color="yellow",
                                        ),
                                    )
                                )
    else: # if callback was triggered by new filtering or no dropdown element was chosen
        fig = get_figure(df, filter_by=xaxis_column_name)
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    fig.update_layout(clickmode='event+select')
    fig.update_traces(marker=dict(size=6,
                                  line=dict(width=1,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers')
                      )
    return fig


@app.callback(
    Output("default-ngl-molecule-dropdown", "value"),
    Input("graph", "clickData"),
)
def update_options(clickData):
    if not clickData:
        raise PreventUpdate

    pt = clickData["points"][0]
    class_label = pt["curveNumber"]
    num = pt["pointNumber"]
    if class_label<len(np.unique(df["class_index"])):
        df_row = df[(df["class_index"] == class_label)].iloc[num]
    else:
        class_label, num_focus = structContainer.get_focus_point()
        df_row = df[(df["class_index"] == class_label)].iloc[num_focus]
    seq_id = df_row['SHORT_ID']
    return seq_id


@app.callback(
    Output("default-ngl-molecule", 'data'),
    Output("default-ngl-molecule", "molStyles"),
    Output("default-ngl-molecule", "stageParameters"),
    Output("default-ngl-molecule", "height"),
    Output("default-ngl-molecule", "width"),
    Input("graph", "clickData"),  # event triggered by click on 3D graph
    # event triggered by dropdown
    Input("default-ngl-molecule-dropdown", "value"),
    Input("nglstyle-dropdown", "value"),
    Input("nglstyle-radio", "value"),
    Input("ngl-stage-quality-dropdown", "value"),
    Input("height-ngl-h", "value"),
    Input("width-ngl-w", "value"),
    Input("chain-atom-input", "n_submit"),
    State('chain-atom-input', 'value'),
    Input("chain-atom-color", "value")
)
def return_molecule(clickData, dropdownVal, style, sidebyside, quality,
                    height, width, res_range_trigger, res_range, color):

    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    if ctx.triggered[0]["prop_id"] == "graph.clickData":
        pt = clickData["points"][0]
        class_label = pt["curveNumber"]
        num = pt["pointNumber"]
        if class_label<len(np.unique(df["class_index"])):
            df_row = df[(df["class_index"] == class_label)].iloc[num]
        else:
            class_label, num_focus = structContainer.get_focus_point()
            df_row = df[(df["class_index"] == class_label)].iloc[num_focus]
        seq_id = df_row['SHORT_ID']
        structContainer.set_structure_ids(seq_id)
    elif ctx.triggered[0]["prop_id"] == "default-ngl-molecule-dropdown.value":
        if len(dropdownVal) == 0:
            raise PreventUpdate
        seq_id = dropdownVal
        structContainer.set_structure_ids(seq_id)
    else:
        pass

    sidebyside_bool = sidebyside == "True"
    molstyles_dict = {
        "representations": style,
        "chosenAtomsColor": color,
        "chosenAtomsRadius": 1,
        "molSpacingXaxis": 100,
        "sideByside": sidebyside_bool
    }
    stage_params = {
        "quality": quality,
    }
    data_p = str(structContainer.get_structur_dir()) + "/"
    seq_ids = structContainer()
    if res_range is not None and (
            bool(re.match(r":[0-9]+-[0-9]+", res_range)) or
            bool(re.match(r":[0-9]+-[0-9]+@[0-9]+", res_range))):
        seq_ids = [seq_id + ".A" + res_range for seq_id in seq_ids]

    data_list = [ngl_parser.get_data(data_path=data_p,  pdb_id=seq_id, color="blue",
                                     reset_view=True, local=True) for seq_id in seq_ids]

    return data_list, molstyles_dict, stage_params, height, width


@app.callback(
    Output("graph-tooltip", "show"),
    Output("graph-tooltip", "bbox"),
    Output("graph-tooltip", "children"),
    Input("graph", "hoverData"),
)
def display_hover(hoverData):
    if hoverData is None:
        return False, no_update, no_update

    pt = hoverData["points"][0]
    bbox = pt["bbox"]
    class_label = pt["curveNumber"]
    num = pt["pointNumber"]
    # if the hovered point is not the focus point
    if class_label<len(np.unique(df["class_index"])):
        df_row = df[(df["class_index"] == class_label)].iloc[num]

    else: # if the hovered point is the focus point (has class_label=n_classes+1)
        print("Focus point")
        class_label, num_focus = structContainer.get_focus_point()
        df_row = df[(df["class_index"] ==class_label)].iloc[num_focus]

    img_src = df_row['IMG_URL']
    #name = df_row['SHORT_ID']

    desc = df_row['hover_data']
    if len(desc) > 300:
        desc = desc[:100] + '...'

    children = [
        html.Div(children=[
            html.Img(src=img_src, style={"width": "100%"}),
            #html.H2(f"{name}", style={"color": "darkblue"}),
            html.P(f"{desc}"),
        ],
            style={'width': '200px', 'whiteSpace': 'normal'})
    ]

    return True, bbox, children


if __name__ == "__main__":
    app.run_server(debug=True, mode='inline')
