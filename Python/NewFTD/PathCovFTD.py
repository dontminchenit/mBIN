import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from scipy.stats import norm
import sys
from numpy import nan
import numpy.ma as ma
import seaborn as sns
from PIL import Image

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from plotNetwork3 import plotNetwork3 
from plotNetwork3_Individual import plotNetwork3Individual
from analysisVisualization import drawCovMatrix, cropImg4
from FTDHelperFunctions import save3DGraph

# Helper Function
def test_corr_sig(corr1, corr2, ncorr1, ncorr2):
    z_corr1 = np.arctanh(corr1)
    z_corr2 = np.arctanh(corr2)
    z_obs = (z_corr1 - z_corr2) / ((1 / (ncorr1 - 3)) + (1 / (ncorr2 - 3))) ** 0.5
    probs = norm.cdf(z_obs)
    return probs

def get_concat_h(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def cmpCoV(cov1, cov2, N, pathNamesRaw, outputDir, cov1_Title, cov1_FileName, cov2_Title, cov2_FileName, cmpCov1_gt_2_Title, cmpCov1_gt_2_FileName, cmpCov2_gt_1_Title, cmpCov2_gt_1_FileName):

    # Covariance matrix for FTD
    covMat1= np.full((N,N), np.nan)
    covMat2= np.full((N,N), np.nan)

    # For Significance difference in Covariance value
    cmpCov1_gt_2 = np.full((N,N), np.nan)
    cmpCov2_gt_1 = np.full((N,N), np.nan)

    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of log %AO for two different anatomical regions for TAU/TDP both not NaN
                N1 = np.sum(~np.isnan(cov1[:,i]) & ~np.isnan(cov1[:,j]))
                N2 = np.sum(~np.isnan(cov2[:,i]) & ~np.isnan(cov2[:,j]))

                if N1 > 3 and N2 >3:
                    # Get the covariance / For only where regions that doesn't have invalid values
                    covMat1[i,j] = ma.corrcoef(ma.masked_invalid(cov1[:,i]), ma.masked_invalid(cov1[:,j]))[0, 1]
                    covMat2[i,j] = ma.corrcoef(ma.masked_invalid(cov2[:,i]), ma.masked_invalid(cov2[:,j]))[0, 1]

                    cmpCov1_gt_2[i,j] = test_corr_sig(covMat2[i,j],covMat1[i,j],N2,N1)<.05 
                    cmpCov2_gt_1[i,j] = test_corr_sig(covMat1[i,j],covMat2[i,j],N1,N2)<.05

    saveName = 'FTD'

    pathNamesAdj = pathNamesRaw.copy() 

    # Clear the figure
    plt.clf()

    # covMat1
    drawCovMatrix(covMat1, pathNamesAdj, pathNamesAdj, cov1_Title, outputDir, cov1_FileName, annot_fontsize = 2)
    # covMat2
    drawCovMatrix(covMat2, pathNamesAdj, pathNamesAdj, cov2_Title, outputDir, cov2_FileName, annot_fontsize = 2)
    # cmpCov1_gt_2
    drawCovMatrix(cmpCov1_gt_2, pathNamesAdj, pathNamesAdj, cmpCov1_gt_2_Title, outputDir, cmpCov1_gt_2_FileName, annot_fontsize = 2)
    # cmpCov2_gt_1
    drawCovMatrix(cmpCov2_gt_1, pathNamesAdj, pathNamesAdj, cmpCov2_gt_1_Title, outputDir, cmpCov2_gt_1_FileName, annot_fontsize = 2)

    return [covMat1, covMat2, cmpCov1_gt_2, cmpCov2_gt_1]


# THERE is no AMY in this dataset
def pathCovFTD(outputDir, NetworkDataGeneral, pathCoM, pathT_GM, pathT_WM, pathNames_3D_Map, sn, plotON = True):

    # ------------------------------------------FTD_3D Atlas Mapping-----------------------------------------
    # Label Names we are able to map to 3D - ['ANG', 'ATC', 'HIP', 'IFC', 'M1', 'MFC', 'OFC', 'PC', 'S1', 'SMTC', 'SPC', 'V1', 'aCING', 'aINS', 'aITC', 'dlPFC', 'iPFC', 'mPFC', 'pCING', 'pSTC']
    # #This MATCHes THE ORDER OF Pathology Regions in the Pathology Dataset (Left to Right)
    LabelNames = pathNames_3D_Map

    # Duplicate LabelNames list with {L, R} Suffix
    suffix_str = ['_L', '_R']
    LabelNames_L = [e + suffix_str[0] for e in LabelNames]
    LabelNames_R = [e + suffix_str[1] for e in LabelNames]

    LabelNames = LabelNames_L + LabelNames_R # Now this matches the order of currCoM too!!
    LabelNames = np.array(LabelNames)

    # Index of Label where (Mean - Marker)Pathology Data is Missing or smaller than prefixed number of Observation --> Therefore is excluded in the map
    TAU_missing_index_GM = []
    TDP_missing_index_GM = []
    TAU_missing_index_WM = []
    TDP_missing_index_WM = []

    # GMWM_list = [pathT_GM, pathT_WM]
    # suffix_GMWM = ['_GM', '_WM']
    
    # ONLY DO GM For NOW
    GMWM_list = [pathT_GM]
    suffix_GMWM = ['_GM']

    for i in range(len(GMWM_list)): # i = 0 : GM / i = 1 : WM
        pathT = GMWM_list[i]
        suffix_M = suffix_GMWM[i]
    
        pathNamesRaw  = list(pathT.columns[5:])

        # Number of Anatomical Regions -> 22 x 2 = 44
        N = len(pathNamesRaw)
        assert(N == 44)

        # Index for the case with tau or tdp for patients with LBD / Pathology Data
        FTD_TAUIndx = (pathT.Tau1_TDP2 == 1)  # False or True
        FTD_TDPIndx = (pathT.Tau1_TDP2 == 2) # False or True

        # Get Log %AO of 22 anatomical regions of the brain
        # aCING, ... --> Name of regions (anatomical region) that does not directly match with Atlas. Therefore we match them later.
        # Masked log for handling the case where the value is NaN
        pathData = np.ma.log(0.01 * pathT.iloc[:, 5:].values + 0.00015).filled(np.nan)

        # Log %AO of FTD TAU vs TDP
        path_TAU = pathData[FTD_TAUIndx,:]
        path_TDP = pathData[FTD_TDPIndx,:]

        # Draw Cov Matrix / CmpCov --> 44 x 44
        covMatlist = cmpCoV(path_TAU, path_TDP, N, pathNamesRaw, outputDir, 'FTD_TAU' + suffix_M,  'CovMat_FTD_TAU' + suffix_M + '.png', 'FTD_TDP' + suffix_M, 'CovMat_FTD_TDP' + suffix_M + '.png', 'FTD - TAU > TDP' + suffix_M, 'FTD_TAU_GT_TDP' + suffix_M + '.png', 'FTD - TDP > TAU' + suffix_M, 'FTD_TDP_GT_TAU' + suffix_M + '.png')

        if plotON:
            # pathCoM and pathToAtlasIndex are ordered by the order of pathNames_3D_Map
            # pathCoM is the Center of mass for 22 anatomical regions A B C / For both L and R Hemisphere
            currCoM = pathCoM

            # Covert (20, 3, 2) --> (40, 3) / first 20 rows {L} and last 20 rows {R} / Order same as pathNames_3D_Map repeated 2 times
            currCoM = np.vstack((currCoM[:, :, 0], currCoM[:, :, 1]))

            # pathNamesRaw - Ordered as same as the Columns in Pathology Dataset (Left to Right)
            # Index of Regions where we can map From the Whole 22 x 2{L, R} = 44 regions
            index_3D_map = [i for i, e in enumerate(pathNamesRaw) if e in LabelNames]

            # Modify the covMatlist using index_3D_map
            for v in range(len(covMatlist)):
                covMatlist[v] = covMatlist[v][np.ix_(index_3D_map, index_3D_map)] # 44 x 44 --> 40 x 40

            # Log %AO of FTD TAU vs TDP
            path_TAU_exp = path_TAU.copy()
            path_TDP_exp = path_TDP.copy()

            # Get only the Part where we can map to 3D - Selecting Specified Columns - sn x 2 Columns
            path_TAU_exp = path_TAU_exp[np.ix_(list(range(path_TAU.shape[0])), index_3D_map)]
            path_TDP_exp = path_TDP_exp[np.ix_(list(range(path_TDP.shape[0])), index_3D_map)]

            # Number of Observation threshold
            obs_thresh = 3

            # Get list containing number of observed subjects per region (list of length 40) / TAU and TDP
            obs_TAU = []
            obs_TDP = []

            for t in range(path_TAU_exp.shape[1]): # iterate over the rows of the 2d array / 40
                non_nan_count = np.count_nonzero(~np.isnan(path_TAU_exp[:, t])) # Number of Non-NaN in this column
                obs_TAU.append(non_nan_count)
            
            for t in range(path_TDP_exp.shape[1]):
                non_nan_count = np.count_nonzero(~np.isnan(path_TDP_exp[:, t])) # Number of Non-NaN in this column
                obs_TDP.append(non_nan_count)

            # TO numpy array
            obs_TAU = np.array(obs_TAU)
            obs_TDP = np.array(obs_TDP)

            print("Label Names")
            print(LabelNames)
            print("number of observations")
            print("TAU")
            print(obs_TAU)
            print("TDP") 
            print(obs_TDP)       
            
            # Save the index of the Pathology Label where the Mean Pathology Value is Missing or less than prefixed number of observations (NaN - TAU, 0 - TDP)
            if i == 0:
                TAU_missing_index_GM = np.argwhere(obs_TAU < obs_thresh).flatten()
                TDP_missing_index_GM = np.argwhere(obs_TDP < obs_thresh).flatten()
            else:
                TAU_missing_index_WM = np.argwhere(obs_TAU < obs_thresh).flatten()
                TDP_missing_index_WM = np.argwhere(obs_TDP < obs_thresh).flatten()

            # Get min/max %AO of LBD
            minPath = np.nanmin(np.vstack([path_TAU_exp, path_TDP_exp]), axis=0)
            maxPath = np.nanmax(np.vstack([path_TAU_exp, path_TDP_exp]) - minPath + 0.0015, axis=0)
            
            # Size of Nodes --> Marker
            MakerVecTAU = np.nanmean(path_TAU_exp, axis=0)
            MakerVecTAU = 3 * (MakerVecTAU - minPath) / maxPath

            MakerVecTDP = np.nanmean(path_TDP_exp, axis=0)
            MakerVecTDP = 3 * (MakerVecTDP - minPath) / maxPath

            cRange = [0, 1]

            # Node color --> Set as red (because cm.jet: 1 --> Red)
            colorVec = np.ones(sn * 2)

            # Images to Crop
            TAU_img = None
            TDP_img = None
            TAU_GT_TDP_img = None
            TDP_GT_TAU_img = None

            covMatNamelist = ['CovMat_FTD_TAU', 'CovMat_FTD_TDP', 'FTD_TAU_GT_TDP', 'FTD_TDP_GT_TAU']
            MakerVecList = [MakerVecTAU, MakerVecTDP]
            imglist = [TAU_img, TDP_img, TAU_GT_TDP_img, TDP_GT_TAU_img]
            for j in range(len(covMatlist)):
                # Define figure
                fig_atlas = plt.figure()

                # Deine MakerVec to Use / j= 0, 2 --> MakerVecYes / j = 1, 3 --> MakerVecNo
                # For cases where we are comparing. We use the first one (ex- TAU gt TDP --> we use TAU MakerVec)
                MakerVec = MakerVecList[j % 2]

                if j == 0 or j == 1:
                    covType = 'original'
                else:
                    covType = 'sig'

                if j == 0 or j == 2: # TAU part
                    currCovMat = np.delete(np.delete(covMatlist[j], TAU_missing_index_GM, axis = 0), TAU_missing_index_GM, axis = 1) # (40 x 40) --> (35 x 35)
                    currPathCoM = np.delete(currCoM, TAU_missing_index_GM, axis = 0) #(40, 3) --> (35, 3)
                    currLabelNames = np.delete(LabelNames, TAU_missing_index_GM, axis=0) # (40,) --> (35,)
                    currMakerVec = np.delete(MakerVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                    currColorVec = np.delete(colorVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                else: # TDP part
                    currCovMat = np.delete(np.delete(covMatlist[j], TDP_missing_index_GM, axis = 0), TDP_missing_index_GM, axis = 1) # (40 x 40) --> (34 x 34)
                    currPathCoM = np.delete(currCoM, TDP_missing_index_GM, axis = 0) #(40, 3) --> (34, 3)
                    currLabelNames = np.delete(LabelNames, TDP_missing_index_GM, axis=0) # (40,) --> (34,)
                    currMakerVec = np.delete(MakerVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)
                    currColorVec = np.delete(colorVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)

                # covMat_TAU, covMat_TDP, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU
                # covMat_TAU[np.ix_(index_3D_map,index_3D_map)] --> Getting parts where we can map to 3D Atlas
                plotNetwork3Individual(fig_atlas, currCovMat, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currPathCoM, currLabelNames, cRange, currMakerVec, currColorVec, 3, showLabels = 0, covType = covType)

                fig_atlas.tight_layout()

                 # Save Graph Object
                save3DGraph(currCovMat, outputDir, covMatNamelist[j] + suffix_M)

                # Save the figure
                plt.savefig(outputDir + '/3D_Atlas_' + covMatNamelist[j] + suffix_M, dpi=1000, bbox_inches='tight')

                # Read figure to crop
                imglist[j] = Image.open(outputDir + '/3D_Atlas_' + covMatNamelist[j] + suffix_M + '.png')
                
                
            comb_img = cropImg4(imglist)
            comb_img.save(outputDir + f'/FTD{suffix_GMWM[i]}.png')

    return LabelNames, TAU_missing_index_GM, TDP_missing_index_GM, TAU_missing_index_WM, TDP_missing_index_WM