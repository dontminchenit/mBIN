import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from plotgiiSurf import plotgiiSurf

def thicknessDist(outputDir, HCResults, PatientTAUResults, PatientTDPResults, NetworkDataGeneral):
    # Assuming you have loaded the data into the following variables:
    # AllResults, graphNames, measureNames, GroupNames, HCIND, svIND, naIND, CoM, LabelLUT, Lausanne250SurfaceMeshFromGII, and numLab

    # Load Data
    i = 4  # thickness
    j = 6  # mean

    # # Measure data
    # measure = AllResults[graphNames[i]][measureNames[j]]
    # if j == 1:  # if betweenness, then we need to normalize
    #     measure = measure / ((numLab - 1) * (numLab - 2))
    
    # Number of (thickness) Regions
    numLab = 400
    
    # Name of Groups we are analyzing
    GroupNames = ['HC', 'TAU', 'TDP']

    # We have the measure for each subgroup (HC, TAU, TDP) - Shape: N x 400
    GroupMeasures = [HCResults['Thickness']['Mean'], PatientTAUResults['Thickness']['Mean'], PatientTDPResults['Thickness']['Mean']]

    # Alpha threshold
    alpha_thresh = 0.05

    # Plotting loop
    testG1All = [1, 2, 1, 2]
    testG2All = [0, 0, 2, 1]

    plotg1 = [0, 0, 3, 0]
    plotg2 = [1, 2, 0, 0]

    for g in range(4): # 1) TAU < HC 2) TDP < HC 3) TAU < TDP 4) TDP < TAU
        G1 = testG1All[g]
        G2 = testG2All[g]
        testName = f"{GroupNames[G1]}<{GroupNames[G2]}"
        TTestsAll_pVal = np.empty(numLab)
        TTestsAll_tstat = np.empty(numLab)

        # Perform t-test on all labels
        for x in range(numLab):
            t, p = ttest_ind(GroupMeasures[G1][:, x], GroupMeasures[G2][:, x], alternative='less') 
            TTestsAll_pVal[x] = p
            TTestsAll_tstat[x] = t

        # Correct for multiple comparison (assuming you have defined the multicmp function or an equivalent) / for alpha value 0.001
        notNaN = ~np.isnan(TTestsAll_pVal) # get only Non-NaN Values of p values
        reject, padj_fdr, alphaSidak, alphaBonf = multipletests(pvals = TTestsAll_pVal[notNaN], alpha = alpha_thresh, method='fdr_bh', is_sorted=False, returnsorted=False)

        # Adjusted p-values/method Name
        padj = padj_fdr
        padjName = 'FDR_Benjamini-Hochberg'

        # allPadjName = ['FDR_Benjamini-Hochberg', 'FWE_Hochberg-Bonferroni', 'FWE_Holm-Bonferroni']

        # Correct p-values for NaN Values
        data = TTestsAll_pVal.copy() # Unadjusted p-values
        data[~notNaN] = 1 # Convert NaN Values to 1
        data[notNaN] = padj # Replace Non-NaN values to adjusted p-values
        
        print(data)
        
        # Center of Mass
        CoM = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['CoM'][0, 0]
        
        # GII Mesh
        GIImesh = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]
        
        # Get only the parts where the p-value is significant
        sig_ind = data < alpha_thresh
        data_sig = data[sig_ind]
        CoM_sig = CoM[sig_ind]

        print(data_sig.shape)
        print(CoM_sig.shape)

        CRange = [0, 1] # crange
        saveName = testName

        # Seperate figure for each Analysis subgroups
        fig = plt.figure()
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

            # Set View Angle
            currAx.view_init(viewE, viewA, viewR)

            # Scatter plot of significant p-values with matching regions
            p = currAx.scatter(CoM_sig[:, 0], CoM_sig[:, 1], CoM_sig[:, 2], c=data_sig, alpha=1, edgecolors='none', cmap='jet', vmin = 0, vmax = alpha_thresh)

            # Plot 3D Brain Mesh
            surfLabels = GIImesh['cdata'][0, 0]
            surfval = GIImesh['cdata'][0, 0]
            surfvalDisp = np.zeros_like(surfval)
            plotgiiSurf(GIImesh['giiSurface_Both'][0,0], surfvalDisp, viewE, viewA, viewR, CRange, 0, 5, currAx, atlasTransparency = 0.01)
        
            # Add colorbar
            fig.colorbar(p, ax=currAx)

            plt.axis('off') # Get rid of axis # Also Get rid of grid and also ticks - Both Rows
            plt.tight_layout() # Makes the mapping tighter --> Bigger

        plt.savefig(os.path.join(outputDir, saveName + '.tif'), dpi=400)


def plotBrain3Points(CoM, GIImesh, data, alphaThresh, panelAll, pSelect_row, cmapTypeForce=None, pointSize=50, surfDisp=None, setViews=None, ViewLoc=None, fillFlag=3):
    
    N = len(data) # 400

    LeftIdx = GIImesh['R0_L1_Index'][0, 0]
    RightIdx = ~GIImesh['R0_L1_Index'][0, 0]

    surfLabels = GIImesh['cdata'][0, 0]
    surfval = GIImesh['cdata'][0, 0]

    XMid = 0
    ZMid = 15

    if surfDisp is not None:
        for i in range(N):
            currLabel = GIImesh['LabelLUT'][i][0]
            currVal = surfDisp[i]
            surfval[surfLabels == currLabel] = currVal

    if alphaThresh > 0:
        PlotTypeRange = 2
        flipCmap = 1
        surfval[surfLabels == 0] = 1
        data[~np.isnan(data)] = np.double(data[~np.isnan(data)] < alphaThresh)
    else:
        PlotTypeRange = 1
        flipCmap = 0

    CoMDisp = CoM.copy()

    if setViews is None or len(setViews) == 0:
        setViews = np.arange(1, 7)
        ViewLoc = np.arange(1, 7)

    if fillFlag is None:
        fillFlag = 3

    for i in range(PlotTypeRange):
        for vSel in range(len(setViews)):  # 6 views L, Lmedial, R, Rmedial, Dorsal, Ventral

            currAx = fig.add_subplot(6, 1, 1, projection='3d')

            v = setViews[vSel]
            if surfDisp is not None:
                valDisp = surfval
            else:
                valDisp = np.zeros(len(surfval))

            regionSel = (data != 0)

            if v == 1:
                viewA = -90
                viewE = 15
                valDisp[RightIdx] = np.nan
                regionSel[CoM[:, 0] > (XMid - 15)] = 0
                lightStr = 2
            elif v == 2:
                viewA = 90
                viewE = 35
                valDisp[RightIdx] = np.nan
                regionSel[(CoM[:, 0] > XMid) | (CoM[:, 0] <= -15)] = 0
                lightStr = 2
            elif v == 3:
                viewA = 90
                viewE = 15
                valDisp[LeftIdx] = np.nan
                regionSel[CoM[:, 0] < (XMid + 15)] = 0
                lightStr = 2
            elif v == 4:
                viewA = -90
                viewE = 35
                valDisp[LeftIdx] = np.nan
                regionSel[(CoM[:, 0] < XMid) | (CoM[:, 0] >= 15)] = 0
                lightStr = 2
            elif v == 5:
                viewA = 0
                viewE = 90
                regionSel[CoM[:, 2] < ZMid] = 0
                lightStr = 1
            elif v == 6:
                viewA = -180
                viewE = -90
                regionSel[CoM[:, 2] >= ZMid] = 0
                lightStr = 1

            if i == 1:  # if normal map
                plotType = None
                cmapType = 1  # jet
                CRange = [0, 1]
            else:  # if threshold map
                valDisp[~np.isnan(valDisp)] = np.double(valDisp[~np.isnan(valDisp)] > alphaThresh)
                plotType = 'thresh'
                CRange = [0, 1]
                cmapType = 2  # gray to red
                flipCmap = 0

            if cmapTypeForce is not None:
                cmapType = cmapTypeForce

            if pointSize is None:
                pointSize = 50

            # if pSelect_row > 0:
            #     panelAll[pSelect_row, ViewLoc[vSel]].select()

            cmap = None
            if cmapType == 1:
                cmap = ListedColormap(plt.cm.Greens(np.linspace(0, 1, 100)))
                cmap.set_under(color='white')
                CRange = [0, 1]
            elif cmapType == 2:
                cmap = ListedColormap(np.flipud(plt.cm.RdBu(np.linspace(0, 1, 100))))
                cmap.set_over(color='white')
                CRange = [-0.25, 0.25]
            elif cmapType == 3:
                cmap = ListedColormap(np.flipud(plt.cm.PiYG(np.linspace(0, 1, 100))))
                cmap.set_over(color='white')
                CRange = [-0.25, 0.25]
            elif cmapType == 4:
                CRange = [0, 6]
                cmap = np.zeros((6, 3))
                cmap[0] = [0, 0, 0]
                cmap[1] = [1, 0, 0]
                cmap[2] = [0.3281, 0.0859, 0.7031]
                cmap[3] = [0, 0.5, 0]
                cmap[4] = [1, 0.5, 0]
                cmap[5] = [1, 1, 0]
                cmap[6] = [240 / 255, 163 / 255, 255 / 255]
            elif cmapType == 5:
                cmap = ListedColormap(plt.cm.Greens(np.linspace(0, 1, 100)))
                cmap.set_over(color='white')
                cmap = cmap(np.arange(cmap.N))
                cmap = np.flipud(cmap)
                flipCmap = 1
            elif cmapType == 6:
                CRange = [-1, 1]
                cmap = np.zeros((5, 3))
                cmap[0] = [45 / 255, 115 / 255, 178 / 255]
                cmap[1] = [106 / 255, 26 / 255, 74 / 255]
                cmap[2] = [1, 1, 1]
                cmap[3] = [0 / 255, 82 / 255, 33 / 255]
                cmap[4] = [228 / 255, 63 / 255, 41 / 255]
            elif cmapType == 7:
                CRange = [0, 6]
                cmap = plt.cm.get_cmap('hot')

            if fillFlag == 1:
                plt.scatter3D(CoMDisp[regionSel, 0], CoMDisp[regionSel, 1], CoMDisp[regionSel, 2], s=pointSize, c=data[regionSel], cmap=cmap, alpha=1, edgecolors='none')
            elif fillFlag == 0:
                plt.scatter3D(CoMDisp[regionSel, 0], CoMDisp[regionSel, 1], CoMDisp[regionSel, 2], s=pointSize, c=data[regionSel], cmap=cmap, alpha=1)
                plt.colorbar()
                t = plotgiiSurf(GIImesh['giiSurface_Both'], valDisp, viewA, viewE, CRange, flipCmap, cmapType, lightStr, cmap)
                t.set_alpha(0.2)
            elif fillFlag == 3:
                
                
                viewR = 0
                currAx.view_init(viewE, viewA, viewR)

                # plt.scatter3D(CoMDisp[regionSel, 0], CoMDisp[regionSel, 1], CoMDisp[regionSel, 2], s=pointSize, c=data[regionSel], cmap=cmap, alpha=1, edgecolors='none')
                # plt.colorbar()
                
                currAx.scatter(CoMDisp[regionSel, 0], CoMDisp[regionSel, 1], CoMDisp[regionSel, 2], c=data[regionSel], cmap=cmap, alpha=1, edgecolors='none')

                plotgiiSurf(GIImesh['giiSurface_Both'][0,0], valDisp, viewE, viewA, viewR, CRange, flipCmap, cmapType, currAx, atlasTransparency = 0.2)
                # def plotgiiSurf(giiSurf, valDisp, viewE, viewA, viewR, cRange, flipCmap, cmapType, currAx, atlasTransparency, cmapIn=None)

                # t = plotgiiSurf(GIImesh['giiSurface_Both'], valDisp, viewA, viewE, CRange, flipCmap, cmapType, lightStr, cmap)
                # t.set_alpha(0.2)
                
            plt.gca().set_visible(False)

# Example usage:
# Replace the input arguments with your actual data.
# plotBrain3Points(CoM, GIImesh, data, alphaThresh, panelAll, pSelect_row, cmapTypeForce, pointSize, surfDisp, setViews, ViewLoc, fillFlag)
# plt.show()


