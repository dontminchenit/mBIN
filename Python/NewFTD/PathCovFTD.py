import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from scipy.stats import norm
import sys
import numpy.ma as ma
import seaborn as sns
from PIL import Image

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from analysisVisualization import drawCovMatrix, nonZeroDegCorr_Path
from FTDHelperFunctions import pathCovGen, pathObsThresh, path3DMapping

def pathCovFTD(outputDir, NetworkDataGeneral, pathCoM, pathT_GM, pathT_WM, plotON = True):
    """Pathology Data Analysis

    Args:
        outputDir (str): Path location to save the analysis results
        NetworkDataGeneral (obj): Preconstructed atlas data
        pathCoM (ndarray): Center of Mass of Pathology Regions (Mean of multiple corresponding Atlas regions' Center of Mass)
        pathT_GM (DataFrame): Pathology %AO data for GM
        pathT_WM (DataFrame): Pathology %AO data for WM
        plotON (bool, optional): Boolean value denoting if we want to do the 3D mapping or not. Defaults to True.

    Returns:
        LabelNames (str list): list of pathology regions that we are able to map to 3D Atlas (currently 20 out of 22 regions)
        TAU_missing_index_GM (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TAU, GM)
        TDP_missing_index_GM (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TDP, GM)
        TAU_missing_index_WM (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TAU, WM)
        TDP_missing_index_WM (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TDP, WM)
    """

    # Label (Pathology Region) Names we are able to map to 3D 
    LabelNames = pathT_GM.columns.values[5:] # Same for both GM and WM

    # Index of Label where (Mean - Marker) Pathology Data is Missing or smaller than prefixed number of Observation --> Therefore is excluded in the 3D mapping
    TAU_missing_index_GM = []
    TDP_missing_index_GM = []
    TAU_missing_index_WM = []
    TDP_missing_index_WM = []

    # GMWM_list = [pathT_GM, pathT_WM] ***
    # suffix_GMWM = ['_GM', '_WM']
    GMWM_list = [pathT_GM] # ONLY DO GM For NOW
    suffix_GMWM = ['_GM']

    for i in range(len(GMWM_list)): # i = 0 : GM / i = 1 : WM
        pthresh_list = [0.05]
        for x in range(len(pthresh_list)):
            # p value threshold
            pthresh = pthresh_list[x]

            # covariance matrix threshold
            cov_thresh = 0.1 # just to make sure there is no Noise
       
            # Pathology Dataframe to Analyze ***
            pathT = GMWM_list[i]

            # Suffix denoting if the analysis is on GM or WM
            suffix_M = suffix_GMWM[i]

            # # List of Pathology Regions in our pathT Data. In alphabetical order (_L first and then _R)
            # pathNamesRaw  = np.array(list(pathT.columns[5:]))

            # Index for the case with tau or tdp for patients *** 
            FTD_TAUIndx = (pathT.Tau1_TDP2 == 1)  # False or True
            FTD_TDPIndx = (pathT.Tau1_TDP2 == 2) # False or True

            # Get Log %AO of 22 anatomical regions of the brain ***
            #pathData = np.ma.log(0.01 * pathT.iloc[:, 5:].values + 0.00015).filled(np.nan) # Masked log for handling the case where the value is NaN
            pathData = np.ma.log(pathT.iloc[:, 5:].values + 0.00015).filled(np.nan)

            # Log %AO of FTD TAU vs TDP ***
            path_TAU = pathData[FTD_TAUIndx,:]
            path_TDP = pathData[FTD_TDPIndx,:]

            # Generate Covariance/Cmp Covariance Matrix for Pathology Data [TAU, TDP, TAU_gt_TDP, TDP_gt_TAU] --> 40 x 40
            covMatlist = pathCovGen(path_TAU, path_TDP, pthresh, cov_thresh)

            # pathCoM is the Center of mass for 20 anatomical regions A B C / For both L and R Hemisphere
            # Covert (20, 3, 2) --> (40, 3) / first 20 rows {L} and last 20 rows {R} / Order same as LabelNames repeated 2 times
            currCoM = np.vstack((pathCoM[:, :, 0], pathCoM[:, :, 1]))

            ################################ Exclude regions where there is few observations ################################
            # This is done after we calculate covariance matrix, to account for the case where TAU_missing_index and TDP_missing_index does not match 
            # Makes TAU_gt_TDP, TDP_gt_TAU hard to compute if we exclude these regions in path_TAU, path_TDP before we compute covariance matrix
            
            # Number of Observation threshold
            obs_thresh = 3

            TAU_missing_index, TDP_missing_index = pathObsThresh(path_TAU, path_TDP, obs_thresh)

            print("TESTING")
            print(TAU_missing_index)
            print(TDP_missing_index)
            
            # Modify the covMatlist data so exclude these regions
            # Rows
            covMatlist[0] = np.delete(covMatlist[0], TAU_missing_index, axis = 0) # TAU Cov Mat
            covMatlist[1] = np.delete(covMatlist[1], TDP_missing_index, axis = 0) # TDP Cov Mat
            covMatlist[2] = np.delete(covMatlist[2], TAU_missing_index, axis = 0) # TAU_gt_TDP Cov Mat
            covMatlist[3] = np.delete(covMatlist[3], TDP_missing_index, axis = 0) # TDP_gt_TAU Cov Mat
            covMatlist[4] = np.delete(covMatlist[4], TAU_missing_index, axis = 0) # covTAU_gt_TDP_raw
            covMatlist[5] = np.delete(covMatlist[5], TDP_missing_index, axis = 0) # covTDP_gt_TAU_raw
            # Columns
            covMatlist[0] = np.delete(covMatlist[0], TAU_missing_index, axis = 1) # TAU Cov Mat
            covMatlist[1] = np.delete(covMatlist[1], TDP_missing_index, axis = 1) # TDP Cov Mat
            covMatlist[2] = np.delete(covMatlist[2], TAU_missing_index, axis = 1) # TAU_gt_TDP Cov Mat
            covMatlist[3] = np.delete(covMatlist[3], TDP_missing_index, axis = 1) # TDP_gt_TAU Cov Mat
            covMatlist[4] = np.delete(covMatlist[4], TAU_missing_index, axis = 1) # covTAU_gt_TDP_raw
            covMatlist[5] = np.delete(covMatlist[5], TDP_missing_index, axis = 1) # covTDP_gt_TAU_raw

            # Also modify path_TAU/path_TDP Data to exclude these regions (columns)
            path_TAU = np.delete(path_TAU, TAU_missing_index, axis = 1)
            path_TDP = np.delete(path_TDP, TDP_missing_index, axis = 1)

            # Generate the LabelNames list for each TAU and TDP
            pathNames_TAU = np.delete(LabelNames, TAU_missing_index)
            pathNames_TDP = np.delete(LabelNames, TDP_missing_index)

            # Generate CoM for each TAU and TDP
            CoM_TAU = np.delete(currCoM, TAU_missing_index, axis = 0)
            CoM_TDP = np.delete(currCoM, TDP_missing_index, axis = 0)

            # Save the index of the Pathology Label where the Mean Pathology Value is Missing or less than prefixed number of observations
            if i == 0: # GM
                TAU_missing_index_GM, TDP_missing_index_GM = TAU_missing_index, TDP_missing_index
            else: # WM
                TAU_missing_index_WM, TDP_missing_index_WM = TAU_missing_index, TDP_missing_index
            ################################################################################################################

            # Draw Cov Matrix / CmpCov --> Now 32 x 32 and Save
            drawCovMatrix(covMatlist[0], pathNames_TAU, pathNames_TAU, 'FTD_TAU' + suffix_M, outputDir, 'CovMat_FTD_TAU' + suffix_M + '.png', annot_fontsize = 2) # TAU covMat
            drawCovMatrix(covMatlist[1], pathNames_TDP, pathNames_TDP, 'FTD_TDP' + suffix_M, outputDir, 'CovMat_FTD_TDP' + suffix_M + '.png', annot_fontsize = 2) # TDP cov Mat
            drawCovMatrix(covMatlist[2], pathNames_TAU, pathNames_TAU, 'FTD: TAU > TDP' + suffix_M, outputDir, 'FTD_TAU_GT_TDP' + suffix_M + '.png', annot_fontsize = 2) # TAU_gt_TDP cov Mat
            drawCovMatrix(covMatlist[3], pathNames_TDP, pathNames_TDP, 'FTD - TDP > TAU' + suffix_M, outputDir, 'FTD_TDP_GT_TAU' + suffix_M + '.png', annot_fontsize = 2) # TDP_gt_TAU cov Mat

            # Draw Nodal Strength vs %AO
            # def nonZeroDegCorr_Path(TAUData, TDPData, covMatTAU, covMatTDP, subplot1_title, subplot2_title, x_label, y_label, outputDir, outputName, linear_regression = False):
            nonZeroDegCorr_Path(path_TAU, path_TDP, covMatlist[0], covMatlist[1], 'FTD_TAU' + suffix_M, 'FTD_TDP' + suffix_M, 'Degree', 'Mean Thickness', outputDir, f"FTD_nonZero_degCorr_Path_{suffix_M}.tif", linear_regression = True)
            plt.clf()

            if plotON: # If we want to plot the 3D Mapping
                # Set fixed density value
                fd_val = 40 # Get top 100

                # Log %AO of FTD TAU vs TDP
                path_TAU_exp = path_TAU.copy()
                path_TDP_exp = path_TDP.copy()

                # Get min/max %AO of LBD
                minPath = np.nanmin(np.vstack([path_TAU_exp, path_TDP_exp]), axis=0)
                maxPath = np.nanmax(np.vstack([path_TAU_exp, path_TDP_exp]) - minPath + 0.0015, axis=0)
                
                # Size of Nodes --> Marker
                MarkerVecTAU = np.nanmean(path_TAU_exp, axis=0)
                MarkerVecTAU = 3 * (MarkerVecTAU - minPath) / maxPath

                MarkerVecTDP = np.nanmean(path_TDP_exp, axis=0)
                MarkerVecTDP = 3 * (MarkerVecTDP - minPath) / maxPath

                cRange = [0, 1]

                # Node color --> Set as red (because cm.jet: 1 --> Red)
                colorVecTAU = np.ones(path_TAU_exp.shape[1])
                colorVecTDP = np.ones(path_TDP_exp.shape[1])

                # 3D Atlas Mapping
                path3DMapping(covMatlist, NetworkDataGeneral, CoM_TAU, pathNames_TAU, MarkerVecTAU, colorVecTAU, CoM_TDP, pathNames_TDP, MarkerVecTDP, colorVecTDP, cRange, outputDir, suffix_M, fd_val)

    return LabelNames, TAU_missing_index_GM, TDP_missing_index_GM, TAU_missing_index_WM, TDP_missing_index_WM