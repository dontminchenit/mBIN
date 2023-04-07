import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from nilearn.plotting import plot_surf_stat_map
import networkx as nx

# plotNetwork3(covMat_yesAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectCol,colorVec,3)
def plotNetwork3(adjMtx, GIImesh, CoM, LabelNames, cRange, MakerVec, panelAll, pSelectCol, colorVec, blue0orng1brwn2jet3, surfDisp=None):

    # What is GII cdata?
    surfLabels = GIImesh['cdata'][0, 0]
    surfval = GIImesh['cdata'][0, 0]

    # What is surfDisp? The number of arguments do not match 
    if surfDisp is not None:
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

    for v in range(2):
        panelAll[v, pSelectCol].select()
        if v == 0:
            viewA = -90
            viewE = 15
        elif v == 1:
            viewA = 0
            viewE = 90

        ColorVRGB = np.zeros((len(colorVec), 3))

        cmap = plt.get_cmap('hot')
        c1 = cmap(round(256 * 1 / 6))
        c2 = cmap(round(256 * 2 / 6))
        c3 = cmap(round(256 * 3 / 6))
        c4 = cmap(round(256 * 4 / 6))
        c5 = cmap(round(256 * 5 / 6))
        c6 = cmap(round(256 * 6 / 6))

        ColorVRGB[colorVec == 0, :] = np.tile(c6, (np.sum(colorVec == 0), 1))
        ColorVRGB[colorVec == 1, :] = np.tile(c1, (np.sum(colorVec == 1), 1))
        ColorVRGB[colorVec == 2, :] = np.tile(c2, (np.sum(colorVec == 2), 1))
        ColorVRGB[colorVec == 3, :] = np.tile(c3, (np.sum(colorVec == 3), 1))
        ColorVRGB[colorVec == 4, :] = np.tile(c4, (np.sum(colorVec == 4), 1))
        ColorVRGB[colorVec == 5, :] = np.tile(c5, (np.sum(colorVec == 5), 1))
        ColorVRGB[colorVec == 6, :] = np.tile(c6, (np.sum(colorVec == 6), 1))
        ColorVRGB[colorVec == 7, :] = np.tile(c6, (np.sum(colorVec == 7), 1))

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
        else:
            h = nx.draw_networkx(G, pos=CoM, node_color=ColorVRGB, node_shape='o', node_size=4.25*MakerVec, width=3)
            h.set_edgecolors(G.edges().values())
            plt.clim(cRange)
            plt.set_cmap('jet')

        showLabels = 1
        for i in range(len(LabelNames)):
            if showLabels:
                nx.draw_networkx_labels(G, {i: CoM[i]}, labels={i: LabelNames[i]}, font_size=8)
            else:
                nx.draw_networkx_labels(G, {i: CoM[i]}, labels={i: ''}, font_size=8)

        plt.axis('equal')
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(GIImesh.giiSurface_Both[0], GIImesh.giiSurface_Both[1], GIImesh.giiSurface_Both[2], facecolors=RGB, alpha=0.05)
        ax.view_init(viewA, viewE)
        plt.show()

        plt.axis('off')
        plt.tight_layout()
        plt.savefig('filename.png', bbox_inches='tight', dpi=300)

