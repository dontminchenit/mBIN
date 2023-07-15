import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from scipy.stats import norm
import sys
from numpy import nan
import numpy.ma as ma
import seaborn as sns
import panel as pn

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from plotNetwork3 import plotNetwork3 
from analysisVisualization import drawCovMatrix

# Helper Function
def test_corr_sig(corr1, corr2, ncorr1, ncorr2):
    z_corr1 = np.arctanh(corr1)
    z_corr2 = np.arctanh(corr2)
    z_obs = (z_corr1 - z_corr2) / ((1 / (ncorr1 - 3)) + (1 / (ncorr2 - 3))) ** 0.5
    probs = norm.cdf(z_obs)
    return probs

def pathCovLBD_noAmy(outputDir, NetworkDataGeneral, pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn):

    # pathCoM is the Center of mass for 6 anatomical regions
    # currCoM gets the Left hemishpere - Core representation (we only pull the left or right - collapsing onto the left hemisphere)
    currCoM = pathCoM[:,:,0]

    # Number of Anatomical Regions -> 6
    N = len(pathNamesRaw)

    # Log %AO of LBD with/without AD
    path_yesAD = pathDataGM[LBD_yesADIndx,:]
    path_noAD = pathDataGM[LBD_noADIndx,:]

    # Covariance matrix for LBD with/without AD
    covMat_yesAD = np.full((N,N), np.nan)
    covMat_noAD = np.full((N,N), np.nan)

    # ??? for LBD with/without AD
    cmpCovYes_gt_No = np.full((N,N), np.nan)
    cmpCovNo_gt_Yes = np.full((N,N), np.nan)

    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of log %AO for two different anatomical regions for with/without AD both not NaN
                NYesAD = np.sum(~np.isnan(path_yesAD[:,i]) & ~np.isnan(path_yesAD[:,j]))
                NNoAD = np.sum(~np.isnan(path_noAD[:,i]) & ~np.isnan(path_noAD[:,j]))

                if NYesAD > 3 and NNoAD >3:
                    # Get the covariance / For only where regions that doesn't have invalid values
                    covMat_yesAD[i,j] = ma.corrcoef(ma.masked_invalid(path_yesAD[:,i]), ma.masked_invalid(path_yesAD[:,j]))[0, 1]
                    covMat_noAD[i,j] = ma.corrcoef(ma.masked_invalid(path_noAD[:,i]), ma.masked_invalid(path_noAD[:,j]))[0, 1]

                    cmpCovYes_gt_No[i,j] = test_corr_sig(covMat_noAD[i,j],covMat_yesAD[i,j],NNoAD,NYesAD)<.05 
                    cmpCovNo_gt_Yes[i,j] = test_corr_sig(covMat_yesAD[i,j],covMat_noAD[i,j],NYesAD,NNoAD)<.05

    # Deleting the Amy - specific to the dataset. / Specify the regions we are interested in at the start. 
    covMat_yesAD = np.delete(covMat_yesAD, 1, axis=0)
    covMat_noAD = np.delete(covMat_noAD, 1, axis=0)
    cmpCovYes_gt_No = np.delete(cmpCovYes_gt_No, 1, axis=0)
    cmpCovNo_gt_Yes = np.delete(cmpCovNo_gt_Yes, 1, axis=0)

    covMat_yesAD = np.delete(covMat_yesAD, 1, axis=1)
    covMat_noAD = np.delete(covMat_noAD, 1, axis=1)
    cmpCovYes_gt_No = np.delete(cmpCovYes_gt_No, 1, axis=1)
    cmpCovNo_gt_Yes = np.delete(cmpCovNo_gt_Yes, 1, axis=1)

    currCoM = np.delete(currCoM, 1, axis=0)

    saveName = 'LBD'

    pathNamesAdj = pathNamesRaw.copy() # 'aCING' 'AMY' 'M1' 'MFC' 'SMTC' 'STRI'
    pathNamesAdj.pop(1) # Delete 'AMY'

    # Clear the figure
    plt.clf()

    #sns.set_theme(style="white")
    
    # ---------------------covMat_yesAD--------------------
    drawCovMatrix(covMat_yesAD, pathNamesAdj, pathNamesAdj, 'LBD-yesAD', outputDir, "/CovMat_yesAD.png")

    # ---------------------covMat_noAD--------------------
    drawCovMatrix(covMat_noAD, pathNamesAdj, pathNamesAdj, 'LBD-noAD', outputDir, "/covMat_noAD.png")

    # ---------------------cmpCovNo_gt_Yes--------------------
    drawCovMatrix(cmpCovNo_gt_Yes, pathNamesAdj, pathNamesAdj, 'noAD > yesAD', outputDir, "/" + saveName + "_1TDPGTTAU.png")

    # ---------------------cmpCovYes_gt_No--------------------
    drawCovMatrix(cmpCovYes_gt_No, pathNamesAdj, pathNamesAdj, 'yesAD > noAD', outputDir, "/" + saveName + "_TAUGTTDP.png")

    # ---------------------LBD_SigPlots--------------------
#---
    # Log %AO of LBD with/without AD
    path_yesAD_exp = path_yesAD.copy()
    path_noAD_exp = path_noAD.copy()

    # Get min/max %AO of LBD
    minPath = np.nanmin(np.vstack([path_yesAD_exp, path_noAD_exp]), axis=0)
    maxPath = np.nanmax(np.vstack([path_yesAD_exp, path_noAD_exp]) - minPath + 0.0015, axis=0)
    
    MakerVecYes = np.nanmean(path_yesAD_exp, axis=0)
    MakerVecYes = 3 * (MakerVecYes - minPath) / maxPath

    MakerVecNo = np.nanmean(path_noAD_exp, axis=0)
    MakerVecNo = 3 * (MakerVecNo - minPath) / maxPath

    MakerVecYes = np.delete(MakerVecYes, 1)
    MakerVecNo = np.delete(MakerVecNo, 1)
#---

    cRange = [0, 1]

    colorVec = np.ones(sn-1)

    LabelNames = pathNamesAdj.copy()

    # Define figure
    fig_atlas = plt.figure()
 
    pSelectCol = 1
    plotNetwork3(fig_atlas, covMat_yesAD, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectCol, colorVec, 3)

    pSelectCol=2
    plotNetwork3(fig_atlas, covMat_noAD, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecNo, pSelectCol, colorVec, 3)

    pSelectCol=3
    plotNetwork3(fig_atlas, cmpCovYes_gt_No, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectCol, colorVec, 3)

    pSelectCol=4
    plotNetwork3(fig_atlas, cmpCovNo_gt_Yes, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecNo, pSelectCol, colorVec, 3)

    fig_atlas.tight_layout()

    # Save the figure
    plt.savefig(outputDir + "/3D_Atlas_Multiple.png", dpi=1000, bbox_inches='tight')


    # ----------TESTING----------
    # fig_testing = plt.figure()

    # pSelectCol = 1
    # plotNetwork3(fig_testing, covMat_yesAD, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectCol, colorVec, 3)

    # fig_testing.tight_layout()

    # plt.savefig(outputDir + "/3D_Atlas_Single.png", dpi=1000, bbox_inches='tight')
    # ----------TESTING----------



# Solid mesh --> Map the dots on the surface, not transparent
###