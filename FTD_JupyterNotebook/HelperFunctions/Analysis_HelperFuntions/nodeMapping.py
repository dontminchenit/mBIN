import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D

import os
import pickle
import sys


def nodeMapping(NetworkDataGeneral, pathCoM, labelNames, markerVec, colorVec, outputDir, title,
                 nodeTransparency = 0.3, atlasTransparency = 0.01, showLabels = 1, surfDisp=None):

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

        # # Create a list of node colors from the ColorVRGB array
        node_color = [tuple(c) for c in ColorVRGB]

        # # Create a list of node sizes from the MakerVec array
        node_size = 4.25 * markerVec

        # Set Axis view
        currAx.view_init(viewE, viewA, viewR)

        # Plot the nodes
        currAx.scatter(xs = pathCoM[:,0], ys = pathCoM[:,1], zs = pathCoM[:,2], s=node_size, c=node_color, alpha=nodeTransparency)

        # Plot the Labels
        if showLabels:
            for i in range(len(labelNames)):
                currAx.text(x=pathCoM[:,0][i], y=pathCoM[:,1][i], z=pathCoM[:,2][i], s=labelNames[i], fontsize = 4)
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

    plt.axis('equal')

    # Save the figure
    plt.savefig(outputDir + '/' + title, dpi=1000, bbox_inches='tight')
    