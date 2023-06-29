import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import compress
import numpy.ma as ma
from scipy.stats import norm
from scipy import stats
import os
import scipy
import sys
from numpy import savetxt
from PIL import Image
from sklearn.linear_model import LinearRegression

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from analysisVisualization import drawScatterplot, nonZeroDegCmp, nonZeroDegCorr, drawCovMatrix, cropImg4, cropImg2, cropImg6, drawthicknessboxplot
from plotNetwork3 import plotNetwork3 
from plotNetwork3_Individual import plotNetwork3Individual 
from FTDHelperFunctions import thickSelectNetwork, thicknessCovMatGenerate, thicknessCovMatGenerate_SigEdgesTSave, sexToBinary, generateZScore, generateWScore, save3DGraph, fixedDensity, thickness3DMapping

def thicknessCovFTD(NetworkDataGeneral, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, outputDir, plotON = True):
    """Thickness Data Analysis

    Args:
        NetworkDataGeneral (obj): Preconstructed atlas data
        pathLUT (DataFrame): Look up table matching Atlas region names to Atlas Labels (Index) in NetworkDataGeneral Object
        HCResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of Healthy Control Subjects
        PatientTAUResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TAU Patients
        PatientTDPResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TDP Patients
        outputDir (str): Path location to save the analysis results
        plotON (bool, optional): Boolean value denoting if we should do 3D Atlas Mapping. Defaults to True.
    """
    #-------------------------------------------------- Select Specific Network we want to do Thickness Analysis on --------------------------------------------------
    # List of network names
    networkNames=['DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All']

    # Denote what NetworkNames to use
    # ii = 5 # 5 --> 'Default'
    ii = 8  # 8 --> Include all networkNames
    
    ########################################################### Variable Summary ######################################################################
    # Label_Name2_list: Names of 400 Regions in the Schaffer Atlas (in order with numlab = 400)
    # currNetwork: List of Boolean denoting if Label_Name2_list contains specified networkNames or not / Boolean list of length 400 / in order of Label_Name2_List
    # saveName: Name of the Network we choose
    # regNames: Thickness Region Names that are in the specific network we choose (in order of the 400 regions)
    # pathNames: Pathology Region Names that matches to regNames (in order of the 400 regions)
    # N: Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
    ###################################################################################################################################################
    currNetwork, saveName, regNames, pathNames, N = thickSelectNetwork(networkNames, ii, pathLUT)

    # Get current Center of Mass for specified current Network / i = 5 --> 91 / i = 8 --> 400
    currCoM = (NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['CoM'][0, 0])[currNetwork]


    #------------------------------------------- From the Selected Thickness Network Divide into HC, TAU, and TDP thickness Data -------------------------------------------
    # PatientTAUResults: thickness mean and volume total values for [Patient (TAU) MRI data IDs (26) x lables (numLab = 400 regions in the sch region if we selected ii = 8)]
    # PatientTDPResults: thickness mean and volume total values for [Patient (TDP) MRI data IDs (30) x lables (numLab = 400 regions in the sch region if we selected ii = 8)]
    # HCResults: thickness mean and volume total values for [Control MRI data IDs (54) x lables (numLab = 400 regions in the sch region if we selected ii = 8)]

    # Get Thickness values from HC MRI Data for only current Network [Control MRI data IDs x N] --> 54 HC x 400 regions
    thickHC = HCResults['Thickness']['Mean'][:,currNetwork]
    # Get Thickness values from Patient TAU MRI Data for only current Network [Control MRI data IDs x N] --> 26 TAU x 400 regions
    thickTAU = PatientTAUResults['Thickness']['Mean'][:,currNetwork]
    # Get Thickness values from Patient TDP MRI Data for only current Network [Control MRI data IDs x N] --> 30 TDP x 400 regions
    thickTDP = PatientTDPResults['Thickness']['Mean'][:,currNetwork]

    # Get NORMALIZED(=Volume/ICV) volume values from HC MRI Data for only current Network [Control MRI data IDs x N] --> 54 HC x 400 regions
    volumeHC = HCResults['Volume']['Normalized'][:,currNetwork]
    # Get NORMALIZED(=Volume/ICV) volume values from Patient TAU MRI Data for only current Network [Control MRI data IDs x N] --> 26 TAU x 400 regions
    volumeTAU = PatientTAUResults['Volume']['Normalized'][:,currNetwork]
    # Get NORMALIZED(=Volume/ICV) volume values from Patient TDP MRI Data for only current Network [Control MRI data IDs x N] --> 30 TDP x 400 regions
    volumeTDP = PatientTDPResults['Volume']['Normalized'][:,currNetwork]

    # Sanity Check
    if ii == 8:
        assert(thickHC.shape == (54, 400))
        assert(thickTAU.shape == (26, 400))
        assert(thickTDP.shape == (30, 400))
        assert(volumeHC.shape == (54, 400))
        assert(volumeTAU.shape == (26, 400))
        assert(volumeTDP.shape == (30, 400))
    else: # ii == 5
        assert(thickHC.shape == (54, 91))
        assert(thickTAU.shape == (26, 91))
        assert(thickTDP.shape == (30, 91))
        assert(volumeHC.shape == (54, 91))
        assert(volumeTAU.shape == (26, 91))
        assert(volumeTDP.shape == (30, 91))

    #------------------------------------------- Calculating Z SCORE -------------------------------------------
    # Get Z score for thickness values
    HC_z, TAU_z, TDP_z = generateZScore(thickHC, thickTAU, thickTDP)

    # Get Z score for VOLUME values
    HC_z_Vol, TAU_z_Vol, TDP_z_Vol = generateZScore(volumeHC, volumeTAU, volumeTDP)

    #------------------------------------------- Calculating W SCORE -------------------------------------------
    # List of 54, 26, and 30 Denoting Age
    ageHC = HCResults['Age'] # (54,)
    ageTAU = PatientTAUResults['Age'] # (26,)
    ageTDP = PatientTDPResults['Age'] # (30,)

    assert(ageHC.shape == (thickHC.shape[0],))
    assert(ageTAU.shape == (thickTAU.shape[0],))
    assert(ageTDP.shape == (thickTDP.shape[0],))
    
    # List of 54, 26, and 30 Denoting Sex
    sexHC = HCResults['Sex'] # (54,)
    sexTAU = PatientTAUResults['Sex'] # (26,)
    sexTDP = PatientTDPResults['Sex'] # (30,)

    # Convert Male --> 0 / Female --> 1
    sexHC = np.array(sexToBinary(sexHC)) #  (54,)
    sexTAU = np.array(sexToBinary(sexTAU)) #  (26,)
    sexTDP = np.array(sexToBinary(sexTDP)) #  (30,)

    assert(sexHC.shape == (thickHC.shape[0],))
    assert(sexTAU.shape == (thickTAU.shape[0],))
    assert(sexTDP.shape == (thickTDP.shape[0],))

    # Stack the Age and Sex features into 2D array
    AgeSexHC = np.squeeze(np.dstack((ageHC, sexHC)))# (54, 2)
    AgeSexTAU = np.squeeze(np.dstack((ageTAU, sexTAU)))  # (26, 2)
    AgeSexTDP = np.squeeze(np.dstack((ageTDP, sexTDP))) # (30, 2)

    assert(AgeSexHC.shape == (thickHC.shape[0], 2))
    assert(AgeSexTAU.shape == (thickTAU.shape[0], 2))
    assert(AgeSexTDP.shape == (thickTDP.shape[0], 2))

    # Get W score for thickness values
    HC_w, TAU_w, TDP_w = generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, thickHC, thickTAU, thickTDP)
    # Get W score for VOLUME values
    HC_w_Vol, TAU_w_Vol, TDP_w_Vol = generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, volumeHC, volumeTAU, volumeTDP)
    
    #------------------------------------------- Generate Box plot of Thickness/Z-Score/W-Score Distributions (Mean across subjects) -------------------------------------------
    drawthicknessboxplot(thickHC, thickTAU, thickTDP, ['HC', 'TAU', 'TDP'], 'Mean Thickness', saveName, outputDir, f"FTD_Distribution_{saveName}_Original.tif")
    drawthicknessboxplot(HC_z, TAU_z, TDP_z, ['HC', 'TAU', 'TDP'], 'Mean Z-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_ZScore.tif")
    drawthicknessboxplot(HC_w, TAU_w, TDP_w, ['HC', 'TAU', 'TDP'], 'Mean W-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_WScore.tif")

    # P value threshold
    # pthresh_list = [0.05, 0.01, 0.005, 0.001]
    pthresh_list = [0.05]
    for x in range(len(pthresh_list)):
        # p value threshold
        pthresh = pthresh_list[x]
        cov_thresh = 0.1 # just to make sure there is no Noise

        #------------------------------------------- COMPUTING CovMat [THICKNESS] -------------------------------------------
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU
        covMatList_Original = thicknessCovMatGenerate_SigEdgesTSave(N, thickHC, thickTAU, thickTDP, pthresh, cov_thresh, outputDir, networkNames, regNames, pathNames, ii)
        covMatList_Z = thicknessCovMatGenerate(N, HC_z, TAU_z, TDP_z, pthresh, cov_thresh) # Z score 
        covMatList_W = thicknessCovMatGenerate(N, HC_w, TAU_w, TDP_w, pthresh, cov_thresh) # W score

         #-------------------------------------------COMPUTING CovMat [VOLUME] -------------------------------------------
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Original_Vol = thicknessCovMatGenerate(N, volumeHC, volumeTAU, volumeTDP, pthresh, cov_thresh) 
        covMatList_Z_Vol = thicknessCovMatGenerate(N, HC_z_Vol, TAU_z_Vol, TDP_z_Vol, pthresh, cov_thresh) # Z score 
        covMatList_W_Vol = thicknessCovMatGenerate(N, HC_w_Vol, TAU_w_Vol, TDP_w_Vol, pthresh, cov_thresh) # W score


        #------------------------------------------- Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient) -------------------------------------------
        ################ THICKNESS ################
        nonZeroDegCmp(covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original.tif") # ORIGINAL VALUE
        plt.clf()
        nonZeroDegCmp(covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score.tif") # Z SCORE
        plt.clf()     
        nonZeroDegCmp(covMatList_W[0], covMatList_W[1], covMatList_W[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score.tif") # W SCORE
        # Clear figure
        plt.clf()
        ################ VOLUME ################
        nonZeroDegCmp(covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original_VOLUME.tif") # ORIGINAL VALUE
        plt.clf()
        nonZeroDegCmp(covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score_VOLUME.tif") # Z SCORE
        plt.clf()
        nonZeroDegCmp(covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score_VOLUME.tif") # W SCORE
        plt.clf()

        #------------------------------------------- Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd  -------------------------------------------
        ################ THICKNESS ################
        nonZeroDegCorr(thickHC, thickTAU, thickTDP, covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original.tif", linear_regression = True)
        plt.clf() # ORIGINAL VALUE
        nonZeroDegCorr(HC_z, TAU_z, TDP_z, covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score.tif", linear_regression = True)
        plt.clf() # Z SCORE
        nonZeroDegCorr(HC_w, TAU_w, TDP_w, covMatList_W[0], covMatList_W[1], covMatList_W[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score.tif", linear_regression = True)
        plt.clf() # W SCORE

        ################ VOLUME ################
        nonZeroDegCorr(volumeHC, volumeTAU, volumeTDP, covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Volume', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original_VOLUME.tif", linear_regression = True)
        plt.clf() # ORIGINAL VALUE
        nonZeroDegCorr(HC_z_Vol, TAU_z_Vol, TDP_z_Vol, covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score_VOLUME.tif", linear_regression = True)
        plt.clf() # Z SCORE
        nonZeroDegCorr(HC_w_Vol, TAU_w_Vol, TDP_w_Vol, covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score_VOLUME.tif", linear_regression = True)
        plt.clf() # W SCORE


        #------------------------------------------- Plot various Covariance Matrix -------------------------------------------
        plotCovMat = False
        if plotCovMat: # Draw various Covaraiance Matrix
            # Covariance Matrix of Healthy Control, TAU, TDP
            drawCovMatrix(covMatList_Original[0], regNames, regNames, 'FTD HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_HC.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[1], regNames, regNames, 'FTD TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[2], regNames, regNames, 'FTD TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP.png", annot_bool = False)

            # Covariance Matrix - TAU, TDP comparing to Healthy Control
            drawCovMatrix(covMatList_Original[3], regNames, regNames, 'FTD TAU > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_HC_pthresh_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[4], regNames, regNames, 'FTD TDP > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_HC_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[5], regNames, regNames, 'FTD TAU < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_lt_HC_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[6], regNames, regNames, 'FTD TDP < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_lt_HC_{pthresh}.png", annot_bool = False)

            # Covariance Matrix  - comparing TAU vs TDP
            drawCovMatrix(covMatList_Original[7], regNames, regNames, 'FTD TAU > TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_TDP_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[8], regNames, regNames, 'FTD TDP > TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_TAU_{pthresh}.png", annot_bool = False)
        
        #-------------------------------------------Plot Covariance Matrix (Network) onto 3D Atlas  -------------------------------------------
        # Only for Original Values (because Z-Score generates the same figure)
        if plotON:   
            # Set Fixed Density Value
            fd_val = 100

            # Set cRange
            cRange = [0, 1]

            # Set labelNames
            LabelNames = regNames

            # Set Node color --> Set as red (because cm.jet: 1 --> Red) Same for all HC/TAU/TDP
            colorVec = np.ones(np.sum(currNetwork))
            
            #------------------------------------------- Original thickness values 3D Mapping -------------------------------------------
            # FOR NODE SIZE (Original Thickness Value)
            thick_HC_exp = thickHC
            thick_TAU_exp = thickTAU
            thick_TDP_exp = thickTDP
            
            # Get the MAX/MIN Thickness values
            minThick_HC = np.nanmin(np.nanmean(thick_HC_exp, axis=0))
            minThick_TAU = np.nanmin(np.nanmean(thick_TAU_exp, axis=0))
            minThick_TDP = np.nanmin(np.nanmean(thick_TDP_exp, axis=0))

            vanishing_val = 0.2
            maxThick_HC = np.nanmax(np.nanmean(thick_HC_exp, axis=0) - minThick_HC) + vanishing_val
            maxThick_TAU = np.nanmax(np.nanmean(thick_TAU_exp, axis=0) - minThick_TAU) + vanishing_val
            maxThick_TDP = np.nanmax(np.nanmean(thick_TDP_exp, axis=0) - minThick_TDP) + vanishing_val
            
            # NODE Size
            MarkerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MarkerVecHC = 3 * (1 - ((MarkerVecHC - minThick_HC) / maxThick_HC))

            MarkerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MarkerVecTAU = 3 * (1 - ((MarkerVecTAU - minThick_TAU) / maxThick_TAU))

            MarkerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MarkerVecTDP = 3 * (1 - ((MarkerVecTDP - minThick_TDP) / maxThick_TDP))

            # Clear figure
            plt.clf()

            # 3D mapping of Original Thickness Values
            thickness3DMapping(covMatList_Original, NetworkDataGeneral, currCoM, LabelNames, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVec, outputDir, pthresh, fd_val, FD = True, W_Score = False)


            #------------------------------------------- W_Score thickness values 3D Mapping -------------------------------------------
            # FOR NODE SIZE (Original Thickness Value)
            thick_HC_exp = HC_w_Vol
            thick_TAU_exp = TAU_w_Vol
            thick_TDP_exp = TDP_w_Vol
            
            # Get the MAX/MIN Thickness values
            minThick_HC = np.nanmin(np.nanmean(thick_HC_exp, axis=0))
            minThick_TAU = np.nanmin(np.nanmean(thick_TAU_exp, axis=0))
            minThick_TDP = np.nanmin(np.nanmean(thick_TDP_exp, axis=0))

            vanishing_val = 0.2
            maxThick_HC = np.nanmax(np.nanmean(thick_HC_exp, axis=0) - minThick_HC) + vanishing_val
            maxThick_TAU = np.nanmax(np.nanmean(thick_TAU_exp, axis=0) - minThick_TAU) + vanishing_val
            maxThick_TDP = np.nanmax(np.nanmean(thick_TDP_exp, axis=0) - minThick_TDP) + vanishing_val
            

            MarkerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MarkerVecHC = 3 * (1 - ((MarkerVecHC - minThick_HC) / maxThick_HC))

            MarkerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MarkerVecTAU = 3 * (1 - ((MarkerVecTAU - minThick_TAU) / maxThick_TAU))

            MarkerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MarkerVecTDP = 3 * (1 - ((MarkerVecTDP - minThick_TDP) / maxThick_TDP))

            # Clear figure
            plt.clf()

            # 3D mapping of W_Score Thickness Values
            thickness3DMapping(covMatList_W_Vol, NetworkDataGeneral, currCoM, LabelNames, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVec, outputDir, pthresh, fd_val, FD = True, W_Score = True)



        