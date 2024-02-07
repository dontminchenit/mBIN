import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

import os
import pickle
import sys

def save3DGraph(adjMtx, outputDir, filename):
    """save Graph as obj that we are mapping to 3D

    Args:
        adjMtx (ndarray): adjacency matrix we use to generate the Graph obj
        outputDir (str): Path location to save the analysis results
        filename (str): Name of the file we are saving
    """
    # Substitute NaN values to 0
    adjMtx[np.isnan(adjMtx)] = 0

    # Get upper triangle of adjMtx
    adjMtx = np.triu(adjMtx)

    # Get graph of adjMtx
    G = nx.Graph(adjMtx)

    # save graph object to file
    pickle.dump(G, open(outputDir + '/' + filename + '.pickle', 'wb'))

    return G
    

# atlasMapping.atlasMapping(NetworkDataGeneral, fd.fixedDensity(cov_volAtPath_w_dict_Drop['TDP_gt_TAU_raw'], fd_val), 
#                           CoM_TDP_Drop, 
#                           pathNames_TDP_Drop, 
#                           markerVecVolTDP, 
#                           colorVecVolTDP, 
#                           thickAtPath_Fig, 'temp',
#                           covType='sig', 
#                           nodeTransparency = 0.3, edgeTransparency = 0.8, atlasTransparency = 0.01, showLabels = 1, surfDisp=None)

def atlasMapping(NetworkDataGeneral, covMat, pathCoM, labelNames, markerVec, colorVec, outputDir, title,
                 covType='original', nodeTransparency = 0.3, edgeTransparency = 0.8, atlasTransparency = 0.01, 
                 showLabels = 1, surfDisp=None):

    # Define figure
    fig = plt.figure()
    
    # Define GII
    GIImesh = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]

    # Define surfLabels and surfval
    surfLabels = GIImesh['cdata'][0, 0]
    surfval = GIImesh['cdata'][0, 0]
    
    # What is surfDisp --> For displaying the surface of the 3D Mesh / Not using for now
    if surfDisp is not None: # NOT USING FOR NOW
        for i in range(len(surfDisp)):
            currLabel = GIImesh['LabelLUT'][i][0] # in MatlabOpaque format --> Cannot parse
            currVal = surfDisp[i]
            surfval[surfLabels == currLabel] = currVal
        surfvalDisp = surfval
    else:
        surfvalDisp = np.zeros_like(surfval)
        
    # Substitute NaN values in covMat to 0
    covMat[np.isnan(covMat)] = 0

    # Get upper triangle of covMat
    covMat = np.triu(covMat)

    # Get graph of covMat
    G = nx.Graph(covMat)
    
    for v in range(2):
        # Get the plt subplot using v(1st, 2nd row) 
        currAx = fig.add_subplot(2, 1, v+1, projection='3d')

        # Set the 3D view angle
        if v == 0:
            viewE = 15
            viewA = 180
            viewR = 0
        elif v == 1:
            viewE = 90
            viewA = -90
            viewR = 0

        ColorVRGB = np.zeros((len(colorVec), 3))

        # Setting the Color Maps - hot
        cmap = plt.get_cmap('hot')

        # Integer index to get the color from a colormap / The values and shape are different from the Matlab code
        c1 = cmap(int(round(256*1/6)))
        c2 = cmap(int(round(256*2/6)))
        c3 = cmap(int(round(256*3/6)))
        c4 = cmap(int(round(256*4/6)))
        c5 = cmap(int(round(256*5/6)))
        c6 = cmap(int(round(256*6/6)))

        # colorVec = np.ones(sn-1)    --> indexed [:3]
        ColorVRGB[colorVec == 0, :] = np.tile(c6[:3], (np.sum(colorVec == 0), 1))
        ColorVRGB[colorVec == 1, :] = np.tile(c1[:3], (np.sum(colorVec == 1), 1)) # --> Only this has any change effect / because: colorVec = np.ones(sn)
        ColorVRGB[colorVec == 2, :] = np.tile(c2[:3], (np.sum(colorVec == 2), 1))
        ColorVRGB[colorVec == 3, :] = np.tile(c3[:3], (np.sum(colorVec == 3), 1))
        ColorVRGB[colorVec == 4, :] = np.tile(c4[:3], (np.sum(colorVec == 4), 1))
        ColorVRGB[colorVec == 5, :] = np.tile(c5[:3], (np.sum(colorVec == 5), 1))
        ColorVRGB[colorVec == 6, :] = np.tile(c6[:3], (np.sum(colorVec == 6), 1))
        ColorVRGB[colorVec == 7, :] = np.tile(c6[:3], (np.sum(colorVec == 7), 1))

        # cmap = 'jet'
        # Create a dictionary of node positions (Center of Mass) from the CoM array
        # {0: array([-9.09673548, 23.98545742, 25.98642635]), 1: array([-43.86865997,  -7.98697188,  41.49483967]), 2: array([-38.51215363,  13.03796983,  37.25460958]), 3: array([-56.05636501, -22.91237593,   0.32856068]), 4: array([-32.29449158,  13.55492001,   2.34797442])}
        pos = {n: pathCoM[i] for i, n in enumerate(G.nodes())}

        # # Create a list of node colors from the ColorVRGB array
        # # [(0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0)]
        node_color = [tuple(c) for c in ColorVRGB]

        # # Create a list of node sizes from the MakerVec array
        # # [7.84399385 7.49494821 6.55088061 8.24452972 7.87500687]
        node_size = 4.25 * markerVec

        # Set the edge colors to the weights of the edges in G
        edge_colors = [G.edges[e]["weight"] for e in G.edges()] # Get weight values of edges in graph G

        # Extract node and edge positions from the layout
        node_xyz = np.array([pos[v] for v in sorted(G)])
        edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])

        currAx.view_init(viewE, viewA, viewR)

        # Plot the nodes
        currAx.scatter(*node_xyz.T, s=node_size, c=node_color, alpha=nodeTransparency)

        # Plot the edges
        for i, vizedge in enumerate(edge_xyz):
            if covType == 'original':
                #colors = plt.cm.jet(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
                # colors = plt.cm.coolwarm(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
                colors = plt.cm.Reds(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
            if covType == 'sig':
#                 colors = 'green'
                colors = 'royalblue'
            currAx.plot(*vizedge.T, c=colors, alpha = edgeTransparency)

        # Plot the Labels
        if showLabels:
            for i in range(len(labelNames)):
                currAx.text(x=pos[i][0], y=pos[i][1], z=pos[i][2], s=labelNames[i], fontsize = 6)
        else:
            pass

        # Plot the 3D Brain surface (Atlas)
        giiSurf = GIImesh['giiSurface_Both'][0,0]
        
        # plots the 3-D triangular surface defined by the points in vectors x, y, and z, and a triangle connectivity matrix T.
        currAx.plot_trisurf(giiSurf['vertices'][0, 0][:,0], giiSurf['vertices'][0, 0][:,1], giiSurf['vertices'][0, 0][:,2],
                            triangles=(giiSurf['faces'][0, 0]) - 1, facecolor='Gray', edgecolor=None, antialiased = False,
                            alpha=atlasTransparency)
        
        # Set Box aspect
        currAx.set_box_aspect([np.ptp(giiSurf['vertices'][0, 0][:, 0]), np.ptp(giiSurf['vertices'][0, 0][:, 1]), 
                               np.ptp(giiSurf['vertices'][0, 0][:, 2])])

        # Set 3D plot angle
        currAx.view_init(viewE, viewA, viewR)
        
        
        
        plt.axis('off') # Get rid of axis # Also Get rid of grid and also ticks - Both Rows
        plt.tight_layout() # Makes the mapping tighter --> Bigger

 
    # plt.axis('off') # Get rid of grid and also ticks - Only get rid of grid for 2nd row
    plt.axis('equal')
    
    # Set title 
    plt.suptitle(title)
    
    # Save the figure
    plt.savefig(outputDir + '/' + title, dpi=1000, bbox_inches='tight')
    
#     # Show Figure
#     plt.show()
    

def graphMapping(NetworkDataGeneral, covMat, pathCoM, labelNames, markerVec, colorVec, outputDir, title,
                 covType='original', nodeTransparency = 0.3, edgeTransparency = 0.8, atlasTransparency = 0.01, 
                 showLabels = 1, surfDisp=None):

    # Define figure
    fig = plt.figure()
    
    # Define GII
    GIImesh = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]

    # Define surfLabels and surfval
    surfLabels = GIImesh['cdata'][0, 0]
    surfval = GIImesh['cdata'][0, 0]
    
    # What is surfDisp --> For displaying the surface of the 3D Mesh / Not using for now
    if surfDisp is not None: # NOT USING FOR NOW
        for i in range(len(surfDisp)):
            currLabel = GIImesh['LabelLUT'][i][0] # in MatlabOpaque format --> Cannot parse
            currVal = surfDisp[i]
            surfval[surfLabels == currLabel] = currVal
        surfvalDisp = surfval
    else:
        surfvalDisp = np.zeros_like(surfval)
        
    # Substitute NaN values in covMat to 0
    covMat[np.isnan(covMat)] = 0

    # Get upper triangle of covMat
    covMat = np.triu(covMat)

    # Get graph of covMat
    G = nx.Graph(covMat)
    
    for v in range(2):
        # Get the plt subplot using v(1st, 2nd row) 
        currAx = fig.add_subplot(2, 1, v+1, projection='3d')

        # Set the 3D view angle
        if v == 0:
            viewE = 15
            viewA = 180
            viewR = 0
        elif v == 1:
            viewE = 90
            viewA = -90
            viewR = 0

        ColorVRGB = np.zeros((len(colorVec), 3))

        # Setting the Color Maps - hot
        cmap = plt.get_cmap('hot')

        # Integer index to get the color from a colormap / The values and shape are different from the Matlab code
        c1 = cmap(int(round(256*1/6)))
        c2 = cmap(int(round(256*2/6)))
        c3 = cmap(int(round(256*3/6)))
        c4 = cmap(int(round(256*4/6)))
        c5 = cmap(int(round(256*5/6)))
        c6 = cmap(int(round(256*6/6)))

        # colorVec = np.ones(sn-1)    --> indexed [:3]
        ColorVRGB[colorVec == 0, :] = np.tile(c6[:3], (np.sum(colorVec == 0), 1))
        ColorVRGB[colorVec == 1, :] = np.tile(c1[:3], (np.sum(colorVec == 1), 1)) # --> Only this has any change effect / because: colorVec = np.ones(sn)
        ColorVRGB[colorVec == 2, :] = np.tile(c2[:3], (np.sum(colorVec == 2), 1))
        ColorVRGB[colorVec == 3, :] = np.tile(c3[:3], (np.sum(colorVec == 3), 1))
        ColorVRGB[colorVec == 4, :] = np.tile(c4[:3], (np.sum(colorVec == 4), 1))
        ColorVRGB[colorVec == 5, :] = np.tile(c5[:3], (np.sum(colorVec == 5), 1))
        ColorVRGB[colorVec == 6, :] = np.tile(c6[:3], (np.sum(colorVec == 6), 1))
        ColorVRGB[colorVec == 7, :] = np.tile(c6[:3], (np.sum(colorVec == 7), 1))

        # cmap = 'jet'
        # Create a dictionary of node positions (Center of Mass) from the CoM array
        # {0: array([-9.09673548, 23.98545742, 25.98642635]), 1: array([-43.86865997,  -7.98697188,  41.49483967]), 2: array([-38.51215363,  13.03796983,  37.25460958]), 3: array([-56.05636501, -22.91237593,   0.32856068]), 4: array([-32.29449158,  13.55492001,   2.34797442])}
        pos = {n: pathCoM[i] for i, n in enumerate(G.nodes())}

        # # Create a list of node colors from the ColorVRGB array
        # # [(0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0)]
        node_color = [tuple(c) for c in ColorVRGB]

        # # Create a list of node sizes from the MakerVec array
        # # [7.84399385 7.49494821 6.55088061 8.24452972 7.87500687]
        node_size = 4.25 * markerVec

        # Set the edge colors to the weights of the edges in G
        edge_colors = [G.edges[e]["weight"] for e in G.edges()] # Get weight values of edges in graph G

        # Extract node and edge positions from the layout
        node_xyz = np.array([pos[v] for v in sorted(G)])
        edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])

        currAx.view_init(viewE, viewA, viewR)

        # Plot the nodes
        currAx.scatter(*node_xyz.T, s=node_size, c=node_color, alpha=nodeTransparency)

        # Plot the edges
        for i, vizedge in enumerate(edge_xyz):
            if covType == 'original':
                #colors = plt.cm.jet(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
                colors = plt.cm.coolwarm(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
            if covType == 'sig':
#                 colors = 'green'
                colors = 'royalblue'
            currAx.plot(*vizedge.T, c=colors, alpha = edgeTransparency)

        # Plot the Labels
        if showLabels:
            for i in range(len(labelNames)):
                currAx.text(x=pos[i][0], y=pos[i][1], z=pos[i][2], s=labelNames[i], fontsize = 6)
        else:
            pass

        # # Plot the 3D Brain surface (Atlas)
        # giiSurf = GIImesh['giiSurface_Both'][0,0]
        
        # # plots the 3-D triangular surface defined by the points in vectors x, y, and z, and a triangle connectivity matrix T.
        # currAx.plot_trisurf(giiSurf['vertices'][0, 0][:,0], giiSurf['vertices'][0, 0][:,1], giiSurf['vertices'][0, 0][:,2],
        #                     triangles=(giiSurf['faces'][0, 0]) - 1, facecolor='Gray', edgecolor=None, antialiased = False,
        #                     alpha=atlasTransparency)
        
        # Set Box aspect
        currAx.set_box_aspect([np.ptp(giiSurf['vertices'][0, 0][:, 0]), np.ptp(giiSurf['vertices'][0, 0][:, 1]), 
                               np.ptp(giiSurf['vertices'][0, 0][:, 2])])

        # Set 3D plot angle
        currAx.view_init(viewE, viewA, viewR)
        
        
        
        plt.axis('off') # Get rid of axis # Also Get rid of grid and also ticks - Both Rows
        plt.tight_layout() # Makes the mapping tighter --> Bigger

 
    # plt.axis('off') # Get rid of grid and also ticks - Only get rid of grid for 2nd row
    plt.axis('equal')
    
    # Set title 
    plt.suptitle(title)
    
    # Save the figure
    plt.savefig(outputDir + '/' + title, dpi=1000, bbox_inches='tight')
    
#     # Show Figure
#     plt.show()