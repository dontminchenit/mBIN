import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from plotgiiSurf import plotgiiSurf

# plotNetwork3(covMat_yesAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectCol,colorVec,3)
def plotNetwork3(currFig, adjMtx, GIImesh, CoM, LabelNames, cRange, MakerVec, pSelectCol, colorVec, blue0orng1brwn2jet3, surfDisp=None):

    # What is GII cdata? [327684 x 1]
    surfLabels = GIImesh['cdata'][0, 0]
    surfval = GIImesh['cdata'][0, 0]

    # What is surfDisp? The number of arguments do not match 
    if surfDisp is not None: # NOT USING FOR NOW
        for i in range(len(surfDisp)):
            currLabel = GIImesh['LabelLUT'][i][0] # in MatlabOpaque format --> Cannot parse
            currVal = surfDisp[i]
            surfval[surfLabels == currLabel] = currVal
        surfvalDisp = surfval
    else:
        surfvalDisp = np.zeros_like(surfval)

    # Substitute NaN values to 0
    adjMtx[np.isnan(adjMtx)] = 0

    # Get upper triangle of adjMtx
    adjMtx = np.triu(adjMtx)

    # Get graph of adjMtx
    G = nx.Graph(adjMtx)

    plotNames = ['L-Lateral', 'Dorsal']

    # Define figure    
    fig = currFig
    
    # Node, Edge Transparency
    nodeTransparency = 0.8
    edgeTransparency = 0.8
    atlasTransparency = 0.01

    for v in range(2):
        # Get the plt subplot using v(1st, 2nd row) pSelectCol / 1: 1, 5    2: 2, 6    3: 3, 7    4: 4, 8
        currAx = fig.add_subplot(2, 4, pSelectCol + (v * 4), projection='3d')

        # Set the 3D view angle
        if v == 0:
            viewE = 15
            viewA = 180
            viewR = 0
        elif v == 1:
            viewE = 90
            viewA = -90
            viewR = 0

        ColorVRGB = np.zeros((len(colorVec), 3)) # np.zeros(5, 3)

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
        ColorVRGB[colorVec == 1, :] = np.tile(c1[:3], (np.sum(colorVec == 1), 1)) # --> Only this has any change effect
        ColorVRGB[colorVec == 2, :] = np.tile(c2[:3], (np.sum(colorVec == 2), 1))
        ColorVRGB[colorVec == 3, :] = np.tile(c3[:3], (np.sum(colorVec == 3), 1))
        ColorVRGB[colorVec == 4, :] = np.tile(c4[:3], (np.sum(colorVec == 4), 1))
        ColorVRGB[colorVec == 5, :] = np.tile(c5[:3], (np.sum(colorVec == 5), 1))
        ColorVRGB[colorVec == 6, :] = np.tile(c6[:3], (np.sum(colorVec == 6), 1))
        ColorVRGB[colorVec == 7, :] = np.tile(c6[:3], (np.sum(colorVec == 7), 1))

        # Draw the network graph and Set the color of the graph edges
        if blue0orng1brwn2jet3 == 1:
            h = nx.draw_networkx(G, pos=CoM, node_color=ColorVRGB, node_shape='o', node_size=4.25*MakerVec, edge_color=(251/255, 97/255, 26/255), width=3)
        elif blue0orng1brwn2jet3 == 2:
            h = nx.draw_networkx(G, pos=CoM, node_color=ColorVRGB, node_shape='o', node_size=4.25*MakerVec, edge_color=(61/255, 43/255, 31/255), width=3)
        elif blue0orng1brwn2jet3 == 0:
            h = nx.draw_networkx(G, pos=CoM, node_color=ColorVRGB, node_shape='o', node_size=4.25*MakerVec, width=3)
        elif blue0orng1brwn2jet3 == 4:
            Cdata = np.transpose(colorVec)
            Cdata[Cdata > 2] = 2
            Cdata[Cdata < -2] = -2
            cmin = -2
            cmax = 2
            cmap = plt.cm.jet
            m = len(cmap)
            index = np.fix((Cdata - cmin) / (cmax - cmin) * m) + 1
            RGB = cmap(index)
            h = nx.draw_networkx(G, pos=CoM, node_color=RGB, node_shape='o', node_size=4.25*MakerVec, width=3)
            h.set_edgecolors(G.edges().values())
            plt.clim(cRange)
            plt.set_cmap('jet')
        else: #blue0orng1brwn2jet3 == 3 --> cmap = 'jet'

            # Create a dictionary of node positions (Center of Mass) from the CoM array
            # {0: array([-9.09673548, 23.98545742, 25.98642635]), 1: array([-43.86865997,  -7.98697188,  41.49483967]), 2: array([-38.51215363,  13.03796983,  37.25460958]), 3: array([-56.05636501, -22.91237593,   0.32856068]), 4: array([-32.29449158,  13.55492001,   2.34797442])}
            pos = {n: CoM[i] for i, n in enumerate(G.nodes())}

            # # Create a list of node colors from the ColorVRGB array
            # # [(0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0), (0.48427829434069847, 0.0, 0.0)]
            node_color = [tuple(c) for c in ColorVRGB]

            # # Create a list of node sizes from the MakerVec array
            # # [7.84399385 7.49494821 6.55088061 8.24452972 7.87500687]
            node_size = 4.25 * MakerVec

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
                colors = plt.cm.jet(edge_colors[i]) # jet color / Already color map is between 0 and 1 so no need for cRange
                currAx.plot(*vizedge.T, c=colors, alpha = edgeTransparency)

            # Plot the Labels
            showLabels = 1
            if showLabels:
                for i in range(len(LabelNames)):
                    currAx.text(x=pos[i][0], y=pos[i][1], z=pos[i][2], s=LabelNames[i], fontsize = 6)
            else:
                pass
            
            # Get rid of axis
            currAx.set_axis_off()

        # Plot the 3D Brain surface (Atlas)
        # def plotgiiSurf(giiSurf, valDisp, viewA, viewE, CRange, flipCmap, cmapType, currAx, cmapIn=None):
        plotgiiSurf(GIImesh['giiSurface_Both'][0,0], surfvalDisp, viewE, viewA, viewR, cRange, 0, 5, currAx, atlasTransparency)

        # plt.axis('off')
        # plt.tight_layout()

        # Maybe enable GPU or in lower resolution.
    
    plt.grid(False)
    plt.axis('off')
    plt.axis('equal')

