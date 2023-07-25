import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.spatial import Delaunay

# plotgiiSurf(GIImesh['giiSurface_Both'][0,0], surfvalDisp, viewE, viewA, viewR, cRange, 0, 5, currAx)
def plotgiiSurf(giiSurf, valDisp, viewE, viewA, viewR, cRange, flipCmap, cmapType, currAx, atlasTransparency, cmapIn=None):
    
    # plots the 3-D triangular surface defined by the points in vectors x, y, and z, and a triangle connectivity matrix T.
    currAx.plot_trisurf(giiSurf['vertices'][0, 0][:,0], giiSurf['vertices'][0, 0][:,1], giiSurf['vertices'][0, 0][:,2], triangles=(giiSurf['faces'][0, 0]) - 1, facecolor='Gray', edgecolor=None, antialiased = False, alpha=atlasTransparency)
    #currAx.plot_trisurf(giiSurf['vertices'][0, 0][:,0], giiSurf['vertices'][0, 0][:,1], giiSurf['vertices'][0, 0][:,2], facecolor='Gray', edgecolor=None, antialiased = False, alpha=atlasTransparency)

    # Lightning (NOT IMPLEMENTED YET)
    # lighting = 'flat'
    # for i in range(lightStr):
    #     a4 = ax.view_init(elev=90)
    #     a4 = ax.view_init(elev=-90)
    
    # Set Box aspect
    currAx.set_box_aspect([np.ptp(giiSurf['vertices'][0, 0][:, 0]), np.ptp(giiSurf['vertices'][0, 0][:, 1]), np.ptp(giiSurf['vertices'][0, 0][:, 2])])
    
    # Set 3D plot angle
    currAx.view_init(viewE, viewA, viewR)
    
    # if cmapIn is None: # NOT USED
    #     if cmapType == 1:
    #         c = plt.cm.RdBu_r(np.linspace(0, 1, 100))
    #         c = np.flipud(c)
    
    #     elif cmapType == 2:
    #         c = plt.cm.Greens(np.linspace(0, 1, 100))
    #     elif cmapType == 3:
    #         c = plt.cm.PuOr(np.linspace(0, 1, 100))
    #         c = np.flipud(c)
    #     elif cmapType == 4:
    #         c = plt.cm.Oranges(np.linspace(0, 1, 100))
    #     elif cmapType == 5:
    #         c = plt.cm.jet(np.linspace(0, 1, 100))
    #         c[0, :] = np.array([1, 1, 1])
    # else:
    #     c = cmapIn

    c = cmapIn
    
    if flipCmap:
        c = np.flipud(c)
        c[0, :] = np.array([1, 1, 1])
