import numpy as np
import os
import sys
import pickle
import scipy
import pandas as pd

import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cm import ScalarMappable

from pycirclize import Circos
from pycirclize.parser import Matrix

import networkx as nx

def condition(source, target, name):
    if name == 'canonical':
        # longer_canon
        cond = ((source == 'L23_2(Cingulate)') & (target == 'L56_1(Cingulate)') | 
                (source == 'L23_3d(Paracingulate)') & (target == 'L56_2(Cingulate)') |
                (source == 'L23_3v(Paracingulate)') & (target == 'L56_2(Cingulate)') |

                (source == 'L23_3d(Paracingulate)') & (target == 'L56_1(Cingulate)') | 
                (source == 'L23_3v(Paracingulate)') & (target == 'L56_1(Cingulate)') | 
                (source == 'L23_4(Rectus)') & (target == 'L56_1(Cingulate)') | 
                (source == 'L23_5(Middle Frontal)') & (target == 'L56_1(Cingulate)') |
                (source == 'L23_4(Rectus)') & (target == 'L56_2(Cingulate)') |
                (source == 'L23_5(Middle Frontal)') & (target == 'L56_2(Cingulate)') |
                (source == 'L23_4(Rectus)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L23_4(Rectus)') & (target == 'L56_3v(Paracingulate)') |
                (source == 'L23_5(Middle Frontal)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L23_5(Middle Frontal)') & (target == 'L56_3v(Paracingulate)'))

    elif name == 'noncanonical':
        # long_noncanon
        cond = ((source == 'L23_1(Cingulate)') & (target == 'L56_2(Cingulate)') | 
                (source == 'L23_2(Cingulate)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L23_2(Cingulate)') & (target == 'L56_3v(Paracingulate)') |

                (source == 'L23_1(Cingulate)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L23_1(Cingulate)') & (target == 'L56_3v(Paracingulate)') | 
                (source == 'L23_1(Cingulate)') & (target == 'L56_4(Rectus)') | 
                (source == 'L23_1(Cingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L23_2(Cingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L23_2(Cingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L23_3d(Paracingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L23_3v(Paracingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L23_3d(Paracingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L23_3v(Paracingulate)') & (target == 'L56_5(Middle Frontal)') | ##
                (source == 'L23_1(Cingulate)') & (target == 'L23_3d(Paracingulate)') | 
                (source == 'L23_1(Cingulate)') & (target == 'L23_3v(Paracingulate)') | 
                (source == 'L23_1(Cingulate)') & (target == 'L23_4(Rectus)') | 
                (source == 'L23_1(Cingulate)') & (target == 'L23_5(Middle Frontal)') | 
                (source == 'L23_2(Cingulate)') & (target == 'L23_4(Rectus)') | 
                (source == 'L23_2(Cingulate)') & (target == 'L23_5(Middle Frontal)') | 
                (source == 'L23_3d(Paracingulate)') & (target == 'L23_4(Rectus)') |
                (source == 'L23_3v(Paracingulate)') & (target == 'L23_4(Rectus)') | 
                (source == 'L23_3d(Paracingulate)') & (target == 'L23_5(Middle Frontal)') |
                (source == 'L23_3v(Paracingulate)') & (target == 'L23_5(Middle Frontal)') | ##
                (source == 'L56_1(Cingulate)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L56_1(Cingulate)') & (target == 'L56_3v(Paracingulate)') | 
                (source == 'L56_1(Cingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L56_1(Cingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L56_2(Cingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L56_2(Cingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L56_3d(Paracingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L56_3v(Paracingulate)') & (target == 'L56_4(Rectus)') |
                (source == 'L56_3d(Paracingulate)') & (target == 'L56_5(Middle Frontal)') |
                (source == 'L56_3v(Paracingulate)') & (target == 'L56_5(Middle Frontal)')
                )

    elif name == 'short':
        # short
        cond = ((source == 'L23_1(Cingulate)') & (target == 'L23_2(Cingulate)') | 
                (source == 'L23_2(Cingulate)') & (target == 'L23_3d(Paracingulate)') |
                (source == 'L23_2(Cingulate)') & (target == 'L23_3v(Paracingulate)') |
                (source == 'L56_1(Cingulate)') & (target == 'L56_2(Cingulate)') |
                (source == 'L56_2(Cingulate)') & (target == 'L56_3d(Paracingulate)') |
                (source == 'L56_2(Cingulate)') & (target == 'L56_3v(Paracingulate)') |
                
                (source == 'L23_1(Cingulate)') & (target == 'L56_1(Cingulate)') | 
                (source == 'L23_2(Cingulate)') & (target == 'L56_2(Cingulate)') | 
                (source == 'L23_3d(Paracingulate)') & (target == 'L56_3d(Paracingulate)') | 
                (source == 'L23_3v(Paracingulate)') & (target == 'L56_3v(Paracingulate)') | 
                (source == 'L23_4(Rectus)') & (target == 'L56_4(Rectus)') | 
                (source == 'L23_5(Middle Frontal)') & (target == 'L56_5(Middle Frontal)'))

    else:
        raise ValueError('Condition Not Defined')
    
    return cond

def chordPlot(covMat, data_label, figTitle, fig_type, markerVec, sub_path = False, sub_path_type = '', ax=None):
    # Substitute NaN values in covMat to 0
    covMat[np.isnan(covMat)] = 0

    # Get upper triangle of covMat
    covMat = np.triu(covMat)

    # Get graph of covMat
    G = nx.Graph(covMat)

    # Creating a mapping from old labels (0-11) to new labels (SMI32_Labels)
    mapping = {old_label: new_label for old_label, new_label in zip(range(12), data_label)}

    # Relabeling the nodes
    G = nx.relabel_nodes(G, mapping)

    # Orderd Positioning the nodes in a circle
    pos = { 'L23_5(Middle Frontal)': np.array([1.00000000e+00, 1.97028653e-08]),
            'L23_4(Rectus)': np.array([0.86602539, 0.50000001]),
            'L23_3v(Paracingulate)': np.array([0.49999998, 0.86602545]),
            'L23_3d(Paracingulate)': np.array([-2.36778241e-08,  1.00000000e+00]),
            'L23_2(Cingulate)': np.array([-0.50000003,  0.86602539]),
            'L23_1(Cingulate)': np.array([-0.86602535,  0.50000007]),
            'L56_1(Cingulate)': np.array([-9.99999960e-01, -6.77199095e-08]),
            'L56_2(Cingulate)': np.array([-0.86602541, -0.49999994]),
            'L56_3d(Paracingulate)': np.array([-0.49999988, -0.86602541]),
            'L56_3v(Paracingulate)': np.array([ 3.19584437e-08, -9.99999960e-01]),
            'L56_4(Rectus)': np.array([ 0.49999992, -0.86602541]),
            'L56_5(Middle Frontal)': np.array([ 0.86602533, -0.50000015])}

    # Define the angle in radians for the clockwise rotation (30 degrees)
    angle = np.radians(15)

    # Perform the rotation for each position
    rotated_pos = {}
    for key, value in pos.items():
        x = value[0]
        y = value[1]
        new_x = x * np.cos(angle) - y * np.sin(angle)
        new_y = x * np.sin(angle) + y * np.cos(angle)
        rotated_pos[key] = np.array([new_x, new_y])

    # Drawing the graph
    fig = plt.figure(figsize=(13, 13))
    edges = G.edges()
    weights = [G[u][v]['weight'] for u,v in edges]

    # Ensure markerVec list is as long as the number of nodes
    if len(markerVec) != len(G.nodes()):
        raise ValueError("node_sizes list must have the same length as the number of nodes")

    if fig_type == 'org':

        # # Normalizing edge weights for color mapping
        # weights_normalized = (weights - np.min(weights)) / (np.max(weights) - np.min(weights))
        # colors = plt.cm.coolwarm(weights_normalized)

        # Normalizing edge weights for color mapping between -1 and 1 / already in range -1 to 1
        weights_normalized = np.array(weights)

        colors = plt.cm.coolwarm((weights_normalized + 1) / 2)  # Rescaling to 0-1 for colormap

        nx.draw_networkx_nodes(G, rotated_pos, node_color='skyblue', node_size=markerVec)
        nx.draw_networkx_labels(G, rotated_pos)
        nx.draw_networkx_edges(G, rotated_pos, edge_color=colors, width=3)

        # Create a ScalarMappable and initialize a colorbar
        sm = ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=-1, vmax=1))
        sm.set_array([])  # You can also set an array of the same shape as data to be represented

        # Add colorbar to the plot
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
        cbar.set_label('Correlation Scale', rotation=0, labelpad=-70)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

        plt.title(figTitle)

    elif fig_type == 'cmp':
        # Use the provided ax for plotting instead of creating a new figure
        plt.sca(ax)

        if not sub_path: # Whole
            nx.draw_networkx_nodes(G, rotated_pos, node_color='skyblue', node_size=markerVec)
            nx.draw_networkx_labels(G, rotated_pos)
            nx.draw_networkx_edges(G, rotated_pos, edge_color='red', width=3)

        else: # Subpaths
            # Iterate through each edge to apply conditions
            for u, v, data in G.edges(data=True):
                edge_weight = data['weight']

                # # Define default edge style
                # edge_color = 'black'
                # edge_width = 2

                # Apply condition - Change color and width for edges with weight above a threshold
                if condition(u, v, sub_path_type): # Condition is met
                    # Draw each edge individually
                    nx.draw_networkx_edges(G, rotated_pos, edgelist=[(u, v)], width=3, edge_color='red')

            nx.draw_networkx_nodes(G, rotated_pos, node_color='skyblue', node_size=markerVec)
            nx.draw_networkx_labels(G, rotated_pos)

        # Create a ScalarMappable and initialize a colorbar
        sm = ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=-1, vmax=1))
        sm.set_array([])  # You can also set an array of the same shape as data to be represented

        # Add colorbar to the plot
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
        cbar.set_label('Significant Difference Scale', rotation=0, labelpad=-70)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

        ax.set_title(figTitle)
    
    #plt.show()
    return fig


def chordPlotRaw(covMat, data_label, figTitle, fig_type, markerVec, sub_path = False, sub_path_type = '', ax=None):
    # Substitute NaN values in covMat to 0
    covMat[np.isnan(covMat)] = 0

    # Get upper triangle of covMat
    covMat = np.triu(covMat)

    # Get graph of covMat
    G = nx.Graph(covMat)

    # Creating a mapping from old labels (0-11) to new labels (SMI32_Labels)
    mapping = {old_label: new_label for old_label, new_label in zip(range(12), data_label)}

    # Relabeling the nodes
    G = nx.relabel_nodes(G, mapping)

    # Orderd Positioning the nodes in a circle
    pos = { 'L23_5(Middle Frontal)': np.array([1.00000000e+00, 1.97028653e-08]),
            'L23_4(Rectus)': np.array([0.86602539, 0.50000001]),
            'L23_3v(Paracingulate)': np.array([0.49999998, 0.86602545]),
            'L23_3d(Paracingulate)': np.array([-2.36778241e-08,  1.00000000e+00]),
            'L23_2(Cingulate)': np.array([-0.50000003,  0.86602539]),
            'L23_1(Cingulate)': np.array([-0.86602535,  0.50000007]),
            'L56_1(Cingulate)': np.array([-9.99999960e-01, -6.77199095e-08]),
            'L56_2(Cingulate)': np.array([-0.86602541, -0.49999994]),
            'L56_3d(Paracingulate)': np.array([-0.49999988, -0.86602541]),
            'L56_3v(Paracingulate)': np.array([ 3.19584437e-08, -9.99999960e-01]),
            'L56_4(Rectus)': np.array([ 0.49999992, -0.86602541]),
            'L56_5(Middle Frontal)': np.array([ 0.86602533, -0.50000015])}

    # Define the angle in radians for the clockwise rotation (30 degrees)
    angle = np.radians(15)

    # Perform the rotation for each position
    rotated_pos = {}
    for key, value in pos.items():
        x = value[0]
        y = value[1]
        new_x = x * np.cos(angle) - y * np.sin(angle)
        new_y = x * np.sin(angle) + y * np.cos(angle)
        rotated_pos[key] = np.array([new_x, new_y])

    # Drawing the graph
    fig = plt.figure(figsize=(13, 13))
    edges = G.edges()
    weights = [G[u][v]['weight'] for u,v in edges]

    # Inverse because p_val
    pval_weights = [1 - 2*x for x in weights]

    # # Normalizing edge weights for color mapping
    # weights_normalized = (pval_weights - np.min(pval_weights)) / (np.max(pval_weights) - np.min(pval_weights))
    # colors = plt.cm.Reds(weights_normalized)

    # Normalizing edge weights for color mapping between -1 and 1 / already in range -1 to 1
    weights_normalized = np.array(pval_weights)

    colors = plt.cm.coolwarm((weights_normalized + 1) / 2)  # Rescaling to 0-1 for colormap
    
    if fig_type == 'cmp':
        # Use the provided ax for plotting instead of creating a new figure
        plt.sca(ax)

        if not sub_path: # Whole
            nx.draw_networkx_nodes(G, rotated_pos, node_color='skyblue', node_size=markerVec)
            nx.draw_networkx_labels(G, rotated_pos)
            nx.draw_networkx_edges(G, rotated_pos, edge_color=colors, width=3)


        else: # Subpaths
            # Get the list of edges from the graph
            edge_list = list(G.edges(data=True))

            # Iterate through each edge to apply conditions
            for idx, (u, v, data) in enumerate(edge_list):

                # Apply condition - Change color and width for edges with weight above a threshold
                if condition(u, v, sub_path_type): # Condition is met
                    # Draw each edge individually
                    nx.draw_networkx_edges(G, rotated_pos, edgelist=[(u, v)], edge_color=colors[idx], width=3)

            nx.draw_networkx_nodes(G, rotated_pos, node_color='skyblue', node_size=markerVec)
            nx.draw_networkx_labels(G, rotated_pos)

         # Create a ScalarMappable and initialize a colorbar
        sm = ScalarMappable(cmap=plt.cm.coolwarm.reversed(), norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])

        # Add colorbar to the plot
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
        cbar.set_label('P-Value Scale', rotation=0, labelpad=-70)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

        ax.set_title(figTitle)
    
    #plt.show()
    return fig
