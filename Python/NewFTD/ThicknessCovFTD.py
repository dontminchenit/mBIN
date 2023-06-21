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
from analysisVisualization import drawScatterplot, nonZeroDegCmp, nonZeroDegCorr, drawCovMatrix, cropImg4, cropImg2, cropImg6
from plotNetwork3 import plotNetwork3 
from plotNetwork3_Individual import plotNetwork3Individual 
from FTDHelperFunctions import thicknessCovMatGenerate, thicknessCovMatGenerate_SigEdgesTSave, sexToBinary, generateZScore, generateWScore, drawthicknessboxplot, save3DGraph, fixedDensity

def thicknessCovFTD(NetworkDataGeneral, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, outputDir, plotON = True):
    # List of network names
    networkNames=['DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All']

    # Denothing what NetworkNames to use
    # 5 --> 'Default'
    #ii = 5
    # 8 --> Include all networkNames
    ii = 8 

    # Add similar way to add on the custom boolean vector
    
    ###################################################################################################################################################
    # Label_Name2_list: Names of 400 Regions in the Schaffer Atlas (in order with numlab = 400)
    # currNetwork: List of Boolean denoting if Label_Name2_list contains specified networkNames or not / Boolean list of length 400 / in order of Label_Name2_List
    # pathLUTSorted: Sort the pathology LUT by Label_ID400x7 1 to 400 / LUT matching 400 thickness regions to PathName
    # pathNames: PathNames that matches to specified networkNames as list (in order of the 400 regions)
    # regNamesRaw: Label_Name2 that contains specified networkNames as list (in order of the 400 regions)
    ###################################################################################################################################################
    
    # Due to mat file issues, loading Label_Name2 from txt file / list of length 400 --> Names of 400 Regions in the Schaffer Atlas (in order)
    # Label_Name2: Schaefer400x7 -> LUT / len() = 400
    with open('/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/NetworkDataGeneral_Schaefer400x7_LUT_Label_Name2.txt') as file:
        Label_Name2_list = [line.rstrip() for line in file]
        file.close()
    
    if(ii == 8):
        currNetwork = np.ones(len(Label_Name2_list), dtype=bool)
    else: # NetworkName --> 'Default' (because --> ii = 5)
        # -1 added because of index difference between python and matlab
        # currNetwork: List of Boolean denoting if Label_Name2_list contains networkNames('Default') or not / Boolean list of length 400
        currNetwork = np.array([networkNames[ii - 1] in label_name for label_name in Label_Name2_list], dtype=bool)

    # saveName --> this specific case: 'Default' / 'All'
    saveName = networkNames[ii - 1] # -1 added because of index difference between python and matlab

    # Sort the pathology LUT by Label_ID400x7
    # pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))
    pathLUTSorted = pathLUT.sort_values(by=['Label_ID400x7']) # Sort by label 1 to 400
    
    # Get the PathNames of the current Network ('Default') inside the pathLUTSorted dataset as list / Out of 400 Regions there are 91 regions that contains 'Default'
    pathNames = pathLUTSorted[currNetwork]['PathName'].tolist()

    # Get the Label_Name2 that contains current Network ('Default') as list
    regNamesRaw = list(compress(Label_Name2_list, currNetwork))

    # Length of regNamesRaw (=91 / 400)
    SN = len(regNamesRaw)

    # Change the regNamesRaw (_ to -)
    regNames = []

    for j in range(SN):
        nameSplit = regNamesRaw[j].replace('_', '-')
        regNames.append(nameSplit)

    # Get current Center of Mass for specified current Network / i = 5 --> 91 / i = 8 --> 400
    currCoM = (NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['CoM'][0, 0])[currNetwork]

    # N = sum(currNetwork) / basically same as SN / 91
    N = SN

    # Sanity Check
    assert(N == 400)

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # PatientTAUResults: thickness mean and volume total values for [Patient (TAU) MRI data IDs (26) x lables (numLab = 400 regions in the sch region)]
    # PatientTDPResults: thickness mean and volume total values for [Patient (TDP) MRI data IDs (30) x lables (numLab = 400 regions in the sch region)]
    # HCResults: thickness mean and volume total values for [Control MRI data IDs (54) x lables (numLab)]

    # Get thickness mean values from HC MRI Data for only current Network [Control MRI data IDs x N] --> 54 HC x 400 regions
    thickHC = HCResults['Thickness']['Mean'][:,currNetwork]
    # Get thickness mean values from Patient TAU MRI Data for only current Network [Control MRI data IDs x N] --> 26 TAU x 400 regions
    thickTAU = PatientTAUResults['Thickness']['Mean'][:,currNetwork]
    # Get thickness mean values from Patient TDP MRI Data for only current Network [Control MRI data IDs x N] --> 30 TDP x 400 regions
    thickTDP = PatientTDPResults['Thickness']['Mean'][:,currNetwork]

    # Get NORMALIZED volume values from HC MRI Data for only current Network [Control MRI data IDs x N] --> 54 HC x 400 regions
    volumeHC = HCResults['Volume']['Normalized'][:,currNetwork]
    # Get volume values from Patient TAU MRI Data for only current Network [Control MRI data IDs x N] --> 26 TAU x 400 regions
    volumeTAU = PatientTAUResults['Volume']['Normalized'][:,currNetwork]
    # Get volume values from Patient TDP MRI Data for only current Network [Control MRI data IDs x N] --> 30 TDP x 400 regions
    volumeTDP = PatientTDPResults['Volume']['Normalized'][:,currNetwork]

    # Sanity Check
    assert(thickHC.shape == (54, 400))
    assert(thickTAU.shape == (26, 400))
    assert(thickTDP.shape == (30, 400))
    assert(volumeHC.shape == (54, 400))
    assert(volumeTAU.shape == (26, 400))
    assert(volumeTDP.shape == (30, 400))

    ################### Z SCORE ###################
    # Get Mean/STD of the thickness values for HC / Mean, STD of HC for each 400 regions / 1 x 400
    hc_Mean_Thick = np.nanmean(thickHC, axis=0)
    hc_SD_Thick = np.nanstd(thickHC, axis=0, ddof=0) # ddof parameter is set to 0, which specifies that the divisor should be N, where N is the number of non-NaN elements in the array

    # Get Z score for thickness values
    HC_z, TAU_z, TDP_z = generateZScore(thickHC, thickTAU, thickTDP, hc_Mean_Thick, hc_SD_Thick)

    # Get Mean/STD of the VOLUME values for HC / Mean, STD of HC for each 400 regions / 1 x 400
    hc_Mean_Vol = np.nanmean(volumeHC, axis=0)
    hc_SD_Vol = np.nanstd(volumeHC, axis=0, ddof=0) 

    # Get Z score for VOLUME values
    HC_z_Vol, TAU_z_Vol, TDP_z_Vol = generateZScore(volumeHC, volumeTAU, volumeTDP, hc_Mean_Vol, hc_SD_Vol)

    ################### W SCORE ###################
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
    
    # Generate Box plot of Thickness/Z-Score/W-Score Distributions (Mean across subjects)
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

        # COMPUTING CovMat [THICKNESS]
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU
        covMatList_Original = thicknessCovMatGenerate_SigEdgesTSave(N, thickHC, thickTAU, thickTDP, pthresh, cov_thresh, outputDir, networkNames, regNames, pathNames, ii)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Z = thicknessCovMatGenerate(N, HC_z, TAU_z, TDP_z, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / W score VERSION!!
        covMatList_W = thicknessCovMatGenerate(N, HC_w, TAU_w, TDP_w, pthresh, cov_thresh)

        # COMPUTING CovMat [VOLUME]
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Original_Vol = thicknessCovMatGenerate(N, volumeHC, volumeTAU, volumeTDP, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Z_Vol = thicknessCovMatGenerate(N, HC_z_Vol, TAU_z_Vol, TDP_z_Vol, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / W score VERSION!!
        covMatList_W_Vol = thicknessCovMatGenerate(N, HC_w_Vol, TAU_w_Vol, TDP_w_Vol, pthresh, cov_thresh)


        # Plotting and Saving figure
        ################ THICKNESS ################
        # 3) Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient) - ORIGINAL VALUE
        nonZeroDegCmp(covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original.tif")
        # Clear figure
        plt.clf()
        # 3) Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient) - Z SCORE
        nonZeroDegCmp(covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score.tif")
        # Clear figure
        plt.clf()
        # 3) Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient) - W SCORE
        nonZeroDegCmp(covMatList_W[0], covMatList_W[1], covMatList_W[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score.tif")
        # Clear figure
        plt.clf()

        ################ VOLUME ################
        nonZeroDegCmp(covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original_VOLUME.tif")
        plt.clf()
        nonZeroDegCmp(covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score_VOLUME.tif")
        plt.clf()
        nonZeroDegCmp(covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score_VOLUME.tif")
        plt.clf()

        ################ THICKNESS ################
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / ORIGINAL VALUE
        nonZeroDegCorr(thickHC, thickTAU, thickTDP, covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original.tif", linear_regression = True)
        # Clear figure
        plt.clf()
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / Z SCORE
        nonZeroDegCorr(HC_z, TAU_z, TDP_z, covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score.tif", linear_regression = True)
        # Clear figure
        plt.clf()
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / W SCORE
        nonZeroDegCorr(HC_w, TAU_w, TDP_w, covMatList_W[0], covMatList_W[1], covMatList_W[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score.tif", linear_regression = True)
        # Clear figure
        plt.clf()

        ################ VOLUME ################
        nonZeroDegCorr(volumeHC, volumeTAU, volumeTDP, covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Volume', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original_VOLUME.tif", linear_regression = True)
        plt.clf()
        nonZeroDegCorr(HC_z_Vol, TAU_z_Vol, TDP_z_Vol, covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score_VOLUME.tif", linear_regression = True)
        plt.clf()
        nonZeroDegCorr(HC_w_Vol, TAU_w_Vol, TDP_w_Vol, covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score_VOLUME.tif", linear_regression = True)
        plt.clf()


        # 5) Plot various Covaraiance Matrix 
        if plotON:

            plotCovMat = False

            if plotCovMat: # Draw various Covaraiance Matrix
                # Covariance Matrix of Control, TAU, TDP
                drawCovMatrix(covMatList_Original[0], regNamesRaw, regNamesRaw, 'FTD HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_HC.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[1], regNamesRaw, regNamesRaw, 'FTD TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[2], regNamesRaw, regNamesRaw, 'FTD TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP.png", annot_bool = False)

                # Covariance Matrix - TAU, TDP comparing to Control
                drawCovMatrix(covMatList_Original[3], regNamesRaw, regNamesRaw, 'FTD TAU > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_HC_pthresh_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[4], regNamesRaw, regNamesRaw, 'FTD TDP > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_HC_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[5], regNamesRaw, regNamesRaw, 'FTD TAU < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_lt_HC_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[6], regNamesRaw, regNamesRaw, 'FTD TDP < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_lt_HC_{pthresh}.png", annot_bool = False)

                # Covariance Matrix  - comparing TAU vs TDP
                drawCovMatrix(covMatList_Original[7], regNamesRaw, regNamesRaw, 'FTD TAU > TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_TDP_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[8], regNamesRaw, regNamesRaw, 'FTD TDP > TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_TAU_{pthresh}.png", annot_bool = False)
            

            # 6) Plot Covariance Matrix (Network) onto 3D Atlas / Only for Original Values (because Z-Score generates the same figure)
            cRange = [0, 1]

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
            
            # MakerVec for HC, TAU and TDP --> FOR NODE SIZE (W_SCORE)
            # Pathology - Bigger the dot --> More disease
            # Thickness - We still want the bigger --> More disease
            # But on thickness the thinner the more disease. 
            ######
            # Change it to the same as Pathology (Normalize 0-1). And then 3 * (1 - normalized thickness value)
            # CHECK if this makes sense in 3D Map. 
            ######
            MakerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MakerVecHC = 3 * (1 - ((MakerVecHC - minThick_HC) / maxThick_HC))

            MakerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MakerVecTAU = 3 * (1 - ((MakerVecTAU - minThick_TAU) / maxThick_TAU))

            MakerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MakerVecTDP = 3 * (1 - ((MakerVecTDP - minThick_TDP) / maxThick_TDP))

            # colorVec for TAU and TDP
            # colorVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            # colorVecTDP = np.nanmean(thick_TDP_exp, axis=0)

            # Node color --> Set as red (because cm.jet: 1 --> Red)
            # colorVecTAU = np.ones(np.sum(currNetwork))
            # colorVecTDP = np.ones(np.sum(currNetwork))
            colorVec = np.ones(np.sum(currNetwork))

            # Originally 0 - but using 3 bc 0 is not yet implemented
            displayType = 3

            LabelNames = regNames

            # Clear figure
            plt.clf()

            # Images to Crop
            HC_img = None
            TAU_img = None
            TDP_img = None
            TAU_gt_HC_img = None
            TDP_gt_HC_img = None
            TAU_lt_HC_img = None
            TDP_lt_HC_img = None
            TAU_gt_TDP_img = None
            TDP_gt_TAU_img = None
            TAU_gt_HC_FD_img = None
            TDP_gt_HC_FD_img = None
            TAU_lt_HC_FD_img = None
            TDP_lt_HC_FD_img = None
            TAU_gt_TDP_FD_img = None
            TDP_gt_TAU_FD_img = None
            imglist = [HC_img, TAU_img, TDP_img, TAU_gt_HC_img, TDP_gt_HC_img, TAU_lt_HC_img, TDP_lt_HC_img, TAU_gt_TDP_img, TDP_gt_TAU_img, TAU_gt_HC_FD_img, TDP_gt_HC_FD_img, TAU_lt_HC_FD_img, TDP_lt_HC_FD_img, TAU_gt_TDP_FD_img, TDP_gt_TAU_FD_img]

            covMatlist = covMatList_Original
            # MODIFY the last 6 covMat in covMatlist for FixedDensity
            fd_val = 100 # Get top 100
            for i in range(9, 15):
                covMatlist[i] = fixedDensity(covMatlist[i], fd_val)

            covMatNamelist = ['covMatHC', 'covMatTAU', 'covMatTDP', 'cmpCovTAU_gt_HC', 'cmpCovTDP_gt_HC', 'cmpCovTAU_lt_HC', 'cmpCovTDP_lt_HC', 'cmpCovTAU_gt_TDP', 'cmpCovTDP_gt_TAU', f'cmpCovTAU_gt_HC_FD_{fd_val}', f'cmpCovTDP_gt_HC_FD_{fd_val}', f'cmpCovTAU_lt_HC_FD_{fd_val}', f'cmpCovTDP_lt_HC_FD_{fd_val}', f'cmpCovTAU_gt_TDP_FD_{fd_val}', f'cmpCovTDP_gt_TAU_FD_{fd_val}']
            MakerVecList = [MakerVecHC, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP]
            
            for j in range(len(covMatlist)):
                # Define figure
                fig_atlas = plt.figure()

                # Define MakerVec to Use / Same np.ones(np.sum(currNetwork))
                MakerVec = MakerVecList[j]

                # Edge color colormap
                if j == 0 or j == 1 or j == 2:
                    covType = 'original'
                else:
                    covType = 'sig'

                # [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]
                plotNetwork3Individual(fig_atlas, covMatlist[j], NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVec, colorVec, 3, showLabels = 0, covType = covType)

                fig_atlas.tight_layout()

                # Save Graph Object
                save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}')

                # Save the figure
                plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}', dpi=1000, bbox_inches='tight')

                # Read figure to crop
                imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}' + '.png')

            thickness_comb_img = cropImg6(imglist[3:9])
            thickness_comb_FD_img = cropImg6(imglist[9:])

            thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(Original).png')
            thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(FixedDensity).png')

            # Images to Crop W SCORE
            HC_img_W = None
            TAU_img_W = None
            TDP_img_W = None
            TAU_gt_HC_img_W = None
            TDP_gt_HC_img_W = None
            TAU_lt_HC_img_W = None
            TDP_lt_HC_img_W = None
            TAU_gt_TDP_img_W = None
            TDP_gt_TAU_img_W = None
            TAU_gt_HC_FD_img_W = None
            TDP_gt_HC_FD_img_W = None
            TAU_lt_HC_FD_img_W = None
            TDP_lt_HC_FD_img_W = None
            TAU_gt_TDP_FD_img_W = None
            TDP_gt_TAU_FD_img_W = None
            
            imglist_W = [HC_img_W, TAU_img_W, TDP_img_W, TAU_gt_HC_img_W, TDP_gt_HC_img_W, TAU_lt_HC_img_W, TDP_lt_HC_img_W, TAU_gt_TDP_img_W, TDP_gt_TAU_img_W, TAU_gt_HC_FD_img_W, TDP_gt_HC_FD_img_W, TAU_lt_HC_FD_img_W, TDP_lt_HC_FD_img_W, TAU_gt_TDP_FD_img_W, TDP_gt_TAU_FD_img_W]

            # FOR NODE SIZE (Original Thickness Value)
            thick_HC_exp = volumeHC
            thick_TAU_exp = volumeTAU
            thick_TDP_exp = volumeTDP
            
            # Get the MAX/MIN Thickness values
            minThick_HC = np.nanmin(np.nanmean(thick_HC_exp, axis=0))
            minThick_TAU = np.nanmin(np.nanmean(thick_TAU_exp, axis=0))
            minThick_TDP = np.nanmin(np.nanmean(thick_TDP_exp, axis=0))

            vanishing_val = 0.2
            maxThick_HC = np.nanmax(np.nanmean(thick_HC_exp, axis=0) - minThick_HC) + vanishing_val
            maxThick_TAU = np.nanmax(np.nanmean(thick_TAU_exp, axis=0) - minThick_TAU) + vanishing_val
            maxThick_TDP = np.nanmax(np.nanmean(thick_TDP_exp, axis=0) - minThick_TDP) + vanishing_val
            
            # MakerVec for HC, TAU and TDP --> FOR NODE SIZE (W_SCORE)
            # Pathology - Bigger the dot --> More disease
            # Thickness - We still want the bigger --> More disease
            # But on thickness the thinner the more disease. 
            ######
            # Change it to the same as Pathology (Normalize 0-1). And then 3 * (1 - normalized thickness value)
            # CHECK if this makes sense in 3D Map. 
            ######
            MakerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MakerVecHC = 3 * (1 - ((MakerVecHC - minThick_HC) / maxThick_HC))

            MakerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MakerVecTAU = 3 * (1 - ((MakerVecTAU - minThick_TAU) / maxThick_TAU))

            MakerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MakerVecTDP = 3 * (1 - ((MakerVecTDP - minThick_TDP) / maxThick_TDP))

            # W-Score is for Volume!!
            covMatlist = covMatList_W_Vol
            # MODIFY the last 6 covMat in covMatlist for FixedDensity
            for i in range(9, 15):
                covMatlist[i] = fixedDensity(covMatlist[i], fd_val)

            covMatNamelist = ['covMatHC_W', 'covMatTAU_W', 'covMatTDP_W', 'cmpCovTAU_gt_HC_W', 'cmpCovTDP_gt_HC_W', 'cmpCovTAU_lt_HC_W', 'cmpCovTDP_lt_HC_W', 'cmpCovTAU_gt_TDP_W', 'cmpCovTDP_gt_TAU_W', f'cmpCovTAU_gt_HC_FD_W_{fd_val}', f'cmpCovTDP_gt_HC_FD_W_{fd_val}', f'cmpCovTAU_lt_HC_FD_W_{fd_val}', f'cmpCovTDP_lt_HC_FD_W_{fd_val}', f'cmpCovTAU_gt_TDP_FD_W_{fd_val}', f'cmpCovTDP_gt_TAU_FD_W_{fd_val}']
            MakerVecList = [MakerVecHC, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP, MakerVecTAU, MakerVecTDP]
            
            for j in range(len(covMatlist)):
                # Define figure
                fig_atlas = plt.figure()

                # Define MakerVec to Use / Same np.ones(np.sum(currNetwork))
                MakerVec = MakerVecList[j]

                # Edge color colormap
                if j == 0 or j == 1 or j == 2:
                    covType = 'original'
                else:
                    covType = 'sig'

                # [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]
                plotNetwork3Individual(fig_atlas, covMatlist[j], NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVec, colorVec, 3, showLabels = 0, covType = covType)

                fig_atlas.tight_layout()

                # Save Graph Object
                save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_WScore')

                # Save the figure
                plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}', dpi=1000, bbox_inches='tight')

                # Read figure to crop
                imglist_W[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}' + '.png')

            thickness_comb_img_W = cropImg6(imglist_W[3:9])
            thickness_comb_FD_img_W = cropImg6(imglist_W[9:])

            thickness_comb_img_W.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(WSCORE_Original).png')
            thickness_comb_FD_img_W.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(WSCORE_FixedDensity).png')

        