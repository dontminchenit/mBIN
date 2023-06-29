import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys
from scipy.stats import norm
from PIL import Image
from sklearn.linear_model import LinearRegression

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from analysisVisualization import drawScatterplot, nonZeroDegCmp, nonZeroDegCorr, drawCovMatrix, cropImg4, cropImg2, cropImg6, drawthicknessboxplot
from plotNetwork3 import plotNetwork3 
from plotNetwork3_Individual import plotNetwork3Individual 
from FTDHelperFunctions import thicknessCovMatGenerate, sexToBinary, generateZScore,generateWScore, save3DGraph, fixedDensity, thicknessAtPathFewObs, thicknessAtPath3DMapping

def thicknessAtPath(outputDir, NetworkDataGeneral, pathCoM, LabelNames_Path, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath, HCResults, PatientTAUResults, PatientTDPResults, TAU_missing_index_GM, TDP_missing_index_GM, plotON = True):
    
    # Covert (20, 3, 2) --> (40, 3) / first 20 rows {L} and last 20 rows {R}
    pathCoM = np.vstack((pathCoM[:, :, 0], pathCoM[:, :, 1]))

    # Plotting and Saving figure Name
    saveName = 'At_Path'

    #------------------------------------------- Loading Thickness values at Pathology Regions -------------------------------------------
    # Convert HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath to appropriate shape: 26 x 20 x 2 --> 26 x 40 / 30 x 20 x 2 --> 30 x 40 / 54 x 20 x 2 --> 54 x 40 / Where 40 is IN ORDER TO pathT Columns {L, R}
    HCthicknessAtPath = np.hstack((HCthicknessAtPath[:, :, 0], HCthicknessAtPath[:, :, 1]))# 54 x 40
    TAUthicknessAtPath = np.hstack((TAUthicknessAtPath[:, :, 0], TAUthicknessAtPath[:, :, 1])) # 26 x 40
    TDPthicknessAtPath = np.hstack((TDPthicknessAtPath[:, :, 0], TDPthicknessAtPath[:, :, 1])) # 30 x 40
    
    #------------------------------------------- Loading Normalized Volume values at Pathology Regions -------------------------------------------
    # Convert HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath to appropriate shape:
    HCnormVolumeAtPath = np.hstack((HCnormVolumeAtPath[:, :, 0], HCnormVolumeAtPath[:, :, 1]))# 54 x 40
    TAUnormVolumeAtPath = np.hstack((TAUnormVolumeAtPath[:, :, 0], TAUnormVolumeAtPath[:, :, 1])) # 26 x 40
    TDPnormVolumeAtPath = np.hstack((TDPnormVolumeAtPath[:, :, 0], TDPnormVolumeAtPath[:, :, 1])) # 30 x 40  
    
    # Sanity Check
    assert(HCthicknessAtPath.shape == (54, 40))
    assert(TAUthicknessAtPath.shape == (26, 40))
    assert(TDPthicknessAtPath.shape == (30, 40))
    assert(HCnormVolumeAtPath.shape == (54, 40))
    assert(TAUnormVolumeAtPath.shape == (26, 40))
    assert(TDPnormVolumeAtPath.shape == (30, 40))

    # N = Number of regions we are analyzing for Thickness at Path
    N = pathCoM.shape[0]

    # Sanity check
    assert(N == 40)
    
    #------------------------------------------- Calculating Z SCORE -------------------------------------------
    # Get Z score for thickness values
    HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z = generateZScore(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath)

    # Get Z score for VOLUME values
    HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol = generateZScore(HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath)
    
    # Sanity Check
    assert(HCthicknessAtPath_z.shape == HCthicknessAtPath.shape)
    assert(TAUthicknessAtPath_z.shape == TAUthicknessAtPath.shape)
    assert(TDPthicknessAtPath_z.shape == TDPthicknessAtPath.shape)
    assert(HCthicknessAtPath_z_Vol.shape == HCthicknessAtPath.shape)
    assert(TAUthicknessAtPath_z_Vol.shape == TAUthicknessAtPath.shape)
    assert(TDPthicknessAtPath_z_Vol.shape == TDPthicknessAtPath.shape)

    #------------------------------------------- Calculating W SCORE -------------------------------------------
    # List of 54, 26, and 30 Denoting Age
    ageHC = HCResults['Age'] # (54,)
    ageTAU = PatientTAUResults['Age'] # (26,)
    ageTDP = PatientTDPResults['Age'] # (30,)

    assert(ageHC.shape == (HCthicknessAtPath.shape[0],))
    assert(ageTAU.shape == (TAUthicknessAtPath.shape[0],))
    assert(ageTDP.shape == (TDPthicknessAtPath.shape[0],))
    
    # List of 54, 26, and 30 Denoting Sex
    sexHC = HCResults['Sex'] # (54,)
    sexTAU = PatientTAUResults['Sex'] # (26,)
    sexTDP = PatientTDPResults['Sex'] # (30,)

    # Convert Male --> 0 / Female --> 1
    sexHC = np.array(sexToBinary(sexHC)) #  (54,)
    sexTAU = np.array(sexToBinary(sexTAU)) #  (26,)
    sexTDP = np.array(sexToBinary(sexTDP)) #  (30,)

    assert(sexHC.shape == (HCthicknessAtPath.shape[0],))
    assert(sexTAU.shape == (TAUthicknessAtPath.shape[0],))
    assert(sexTDP.shape == (TDPthicknessAtPath.shape[0],))

    # Stack the Age and Sex features into 2D array
    AgeSexHC = np.squeeze(np.dstack((ageHC, sexHC)))# (54, 2)
    AgeSexTAU = np.squeeze(np.dstack((ageTAU, sexTAU)))  # (26, 2)
    AgeSexTDP = np.squeeze(np.dstack((ageTDP, sexTDP))) # (30, 2)

    assert(AgeSexHC.shape == (HCthicknessAtPath.shape[0], 2))
    assert(AgeSexTAU.shape == (TAUthicknessAtPath.shape[0], 2))
    assert(AgeSexTDP.shape == (TDPthicknessAtPath.shape[0], 2))

    # Get W score for thickness values
    HC_w, TAU_w, TDP_w = generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath)
    # Get W score for VOLUME values
    HC_w_Vol, TAU_w_Vol, TDP_w_Vol = generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath)

    # Sanity Check
    assert(HC_w.shape == HCthicknessAtPath.shape)
    assert(TAU_w.shape == TAUthicknessAtPath.shape)
    assert(TDP_w.shape == TDPthicknessAtPath.shape)
    assert(HC_w_Vol.shape == HCthicknessAtPath.shape)
    assert(TAU_w_Vol.shape == TAUthicknessAtPath.shape)
    assert(TDP_w_Vol.shape == TDPthicknessAtPath.shape)

    #------------------------------------------- Generate Box plot of Thickness/Z-Score/W-Score Distributions (Mean across subjects) -------------------------------------------
    drawthicknessboxplot(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, ['HC', 'TAU', 'TDP'], 'Mean Thickness', saveName, outputDir, f"FTD_Distribution_{saveName}_Original.tif")
    drawthicknessboxplot(HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, ['HC', 'TAU', 'TDP'], 'Mean Z-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_ZScore.tif")
    drawthicknessboxplot(HC_w, TAU_w, TDP_w, ['HC', 'TAU', 'TDP'], 'Mean W-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_WScore.tif")

    # P value threshold
    # pthresh_list = [0.05, 0.01, 0.005, 0.001]
    pthresh_list = [0.05]
    for x in range(len(pthresh_list)):
        # p value threshold
        pthresh = pthresh_list[x]
        cov_thresh = 0.1 # just to make sure there is no Noise
        
        #------------------------------------------- COMPUTING CovMat [THICKNESS] -------------------------------------------
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU, cmpCovTAU_gt_HC_raw, cmpCovTDP_gt_HC_raw, cmpCovTAU_lt_HC_raw, cmpCovTDP_lt_HC_raw, cmpCovTAU_gt_TDP_raw, cmpCovTDP_gt_TAU_raw
        covMatList_Original = thicknessCovMatGenerate(N, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, pthresh, cov_thresh)
        covMatList_Z = thicknessCovMatGenerate(N, HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, pthresh, cov_thresh) # Z score 
        covMatList_W = thicknessCovMatGenerate(N, HC_w, TAU_w, TDP_w, pthresh, cov_thresh)                                              # W score
 
        #------------------------------------------- COMPUTING CovMat [VOLUME] -------------------------------------------
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU, cmpCovTAU_gt_HC_raw, cmpCovTDP_gt_HC_raw, cmpCovTAU_lt_HC_raw, cmpCovTDP_lt_HC_raw, cmpCovTAU_gt_TDP_raw, cmpCovTDP_gt_TAU_raw
        covMatList_Original_Vol = thicknessCovMatGenerate(N, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath, pthresh, cov_thresh)
        covMatList_Z_Vol = thicknessCovMatGenerate(N, HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol, pthresh, cov_thresh) # Z score
        covMatList_W_Vol = thicknessCovMatGenerate(N, HC_w_Vol, TAU_w_Vol, TDP_w_Vol, pthresh, cov_thresh)                                              # W score


        #------------------------------------------- EXCLUDING REGIONS WHERE THERE IS FEW PATHOLOGY OBSERVATIONS -------------------------------------------
        ####### Z Score #######
        TAUthicknessAtPath_z = np.delete(TAUthicknessAtPath_z, TAU_missing_index_GM, axis = 1)
        TDPthicknessAtPath_z = np.delete(TDPthicknessAtPath_z, TDP_missing_index_GM, axis = 1)

        TAUthicknessAtPath_z_Vol = np.delete(TAUthicknessAtPath_z_Vol, TAU_missing_index_GM, axis = 1)
        TDPthicknessAtPath_z_Vol = np.delete(TDPthicknessAtPath_z_Vol, TDP_missing_index_GM, axis = 1)
        
        ####### W Score #######
        TAU_w = np.delete(TAU_w, TAU_missing_index_GM, axis = 1)
        TDP_w = np.delete(TDP_w, TDP_missing_index_GM, axis = 1)

        TAU_w_Vol = np.delete(TAU_w_Vol, TAU_missing_index_GM, axis = 1)
        TDP_w_Vol = np.delete(TDP_w_Vol, TDP_missing_index_GM, axis = 1)

        ####### Covariance Matrix #######
        covMatList_Original = thicknessAtPathFewObs(covMatList_Original, TAU_missing_index_GM, TDP_missing_index_GM)
        covMatList_Z = thicknessAtPathFewObs(covMatList_Z, TAU_missing_index_GM, TDP_missing_index_GM)
        covMatList_W = thicknessAtPathFewObs(covMatList_W, TAU_missing_index_GM, TDP_missing_index_GM)

        covMatList_Original_Vol = thicknessAtPathFewObs(covMatList_Original_Vol, TAU_missing_index_GM, TDP_missing_index_GM)
        covMatList_Z_Vol = thicknessAtPathFewObs(covMatList_Z_Vol, TAU_missing_index_GM, TDP_missing_index_GM)
        covMatList_W_Vol = thicknessAtPathFewObs(covMatList_W_Vol, TAU_missing_index_GM, TDP_missing_index_GM)
        
        ####### Thickness Values #######
        TAUthicknessAtPath = np.delete(TAUthicknessAtPath, TAU_missing_index_GM, axis = 1)
        TDPthicknessAtPath = np.delete(TDPthicknessAtPath, TDP_missing_index_GM, axis = 1)

        ####### Normalized Volume Values #######
        TAUnormVolumeAtPath = np.delete(TAUnormVolumeAtPath, TAU_missing_index_GM, axis = 1)
        TDPnormVolumeAtPath = np.delete(TDPnormVolumeAtPath, TDP_missing_index_GM, axis = 1)

        ######## Generate the LabelNames_Path list for each TAU and TDP #######
        LabelNames_Path_HC = LabelNames_Path
        LabelNames_Path_TAU = np.delete(LabelNames_Path, TAU_missing_index_GM)
        LabelNames_Path_TDP = np.delete(LabelNames_Path, TDP_missing_index_GM)

        ######## Generate CoM for each TAU and TDP ########
        CoM_HC = pathCoM
        CoM_TAU = np.delete(pathCoM, TAU_missing_index_GM, axis = 0)
        CoM_TDP = np.delete(pathCoM, TDP_missing_index_GM, axis = 0)
        
        #------------------------------------------- Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient) - -------------------------------------------
        ################ THICKNESS ################  
        nonZeroDegCmp(covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original.tif") 
        plt.clf() # ORIGINAL VALUE
        nonZeroDegCmp(covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score.tif")
        plt.clf() # Z score
        nonZeroDegCmp(covMatList_W[0], covMatList_W[1], covMatList_W[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Thickness", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score.tif")
        plt.clf() # W score

        ################ VOLUME ################
        nonZeroDegCmp(covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], ['HC', 'TAU', 'TDP'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Original_VOLUME.tif")
        plt.clf()
        nonZeroDegCmp(covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], ['HC_Z', 'TAU_Z', 'TDP_Z'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_Z_Score_VOLUME.tif")
        plt.clf()
        nonZeroDegCmp(covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], ['HC_W', 'TAU_W', 'TDP_W'], 'Degree', saveName + " Volume", outputDir, f"FTD_nonZero_degCmp_{saveName}_W_Score_VOLUME.tif")
        plt.clf()

        #------------------------------------------- Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd -------------------------------------------
        ################ THICKNESS ################
        nonZeroDegCorr(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original.tif", linear_regression = True)
        plt.clf() # Original Value
        nonZeroDegCorr(HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score.tif", linear_regression = True)
        plt.clf() # Z Score
        nonZeroDegCorr(HC_w, TAU_w, TDP_w, covMatList_W[0], covMatList_W[1], covMatList_W[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score.tif", linear_regression = True)
        plt.clf() # W Score

        ################ VOLUME ################
        nonZeroDegCorr(HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath, covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Volume', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original_VOLUME.tif", linear_regression = True)
        plt.clf() # Original Value
        nonZeroDegCorr(HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol, covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score_VOLUME.tif", linear_regression = True)
        plt.clf() # Z Score
        nonZeroDegCorr(HC_w_Vol, TAU_w_Vol, TDP_w_Vol, covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score_VOLUME.tif", linear_regression = True)
        plt.clf() # W Score

        #------------------------------------------- Plot various Covariance Matrix -------------------------------------------
        plotCovMat = False
        if plotCovMat: # Draw various Covaraiance Matrix
            # Covariance Matrix of Control, TAU, TDP
            drawCovMatrix(covMatList_Original[0], LabelNames_Path_HC, LabelNames_Path_HC, 'FTD HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_HC.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[1], LabelNames_Path_TAU, LabelNames_Path_TAU, 'FTD TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[2], LabelNames_Path_TDP, LabelNames_Path_TDP, 'FTD TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP.png", annot_bool = False)

            # Covariance Matrix - TAU, TDP comparing to Control
            drawCovMatrix(covMatList_Original[3], LabelNames_Path_TAU, LabelNames_Path_TAU, 'FTD TAU > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_HC_pthresh_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[4], LabelNames_Path_TDP, LabelNames_Path_TDP, 'FTD TDP > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_HC_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[5], LabelNames_Path_TAU, LabelNames_Path_TAU, 'FTD TAU < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_lt_HC_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[6], LabelNames_Path_TDP, LabelNames_Path_TDP, 'FTD TDP < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_lt_HC_{pthresh}.png", annot_bool = False)

            # Covariance Matrix  - comparing TAU vs TDP
            drawCovMatrix(covMatList_Original[7], LabelNames_Path_TAU, LabelNames_Path_TAU, 'FTD TAU > TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_TDP_{pthresh}.png", annot_bool = False)
            drawCovMatrix(covMatList_Original[8], LabelNames_Path_TDP, LabelNames_Path_TDP, 'FTD TDP > TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_TAU_{pthresh}.png", annot_bool = False)
        
        
        #-------------------------------------------Plot Covariance Matrix (Network) onto 3D Atlas  -------------------------------------------
        # For Original Values (because Z-Score generates the same figure)
        if plotON:
            # Set cRange
            cRange = [0, 1]
            
            # Set fixed density value
            fd_val = 100 # Get top 100
            
            # Node color --> Set as red (because cm.jet: 1 --> Red)
            colorVecHC = np.ones(HCthicknessAtPath.shape[1])
            colorVecTAU = np.ones(TAUthicknessAtPath.shape[1])
            colorVecTDP = np.ones(TDPthicknessAtPath.shape[1])
            
            #------------------------------------------- Original thickness values 3D Mapping -------------------------------------------
            # FOR NODE SIZE (Original Thickness Value)
            thick_HC_exp = HCthicknessAtPath
            thick_TAU_exp = TAUthicknessAtPath
            thick_TDP_exp = TDPthicknessAtPath

            # Get the MAX/MIN Thickness values
            minThick_HC = np.nanmin(np.nanmean(thick_HC_exp, axis=0))
            minThick_TAU = np.nanmin(np.nanmean(thick_TAU_exp, axis=0))
            minThick_TDP = np.nanmin(np.nanmean(thick_TDP_exp, axis=0))

            vanishing_val = 0.2
            maxThick_HC = np.nanmax(np.nanmean(thick_HC_exp, axis=0) - minThick_HC) + vanishing_val
            maxThick_TAU = np.nanmax(np.nanmean(thick_TAU_exp, axis=0) - minThick_TAU) + vanishing_val
            maxThick_TDP = np.nanmax(np.nanmean(thick_TDP_exp, axis=0) - minThick_TDP) + vanishing_val
            
            # MakerVec for HC, TAU and TDP -> FOR NODE SIZE
            MarkerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MarkerVecHC = 3 * (1 - ((MarkerVecHC - minThick_HC) / maxThick_HC))

            MarkerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MarkerVecTAU = 3 * (1 - ((MarkerVecTAU - minThick_TAU) / maxThick_TAU))

            MarkerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MarkerVecTDP = 3 * (1 - ((MarkerVecTDP - minThick_TDP) / maxThick_TDP))

            # Clear figure
            plt.clf()

            thicknessAtPath3DMapping(covMatList_Original, NetworkDataGeneral, CoM_HC, CoM_TAU, CoM_TDP, LabelNames_Path_HC, LabelNames_Path_TAU, LabelNames_Path_TDP, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVecHC, colorVecTAU, colorVecTDP, outputDir, pthresh, fd_val, FD = True, W_Score = False)

            #------------------------------------------- W_Score thickness values 3D Mapping -------------------------------------------
            # FOR NODE SIZE (Normalized Volume)
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

            
            # MarkerVec for HC, TAU and TDP -> FOR NODE SIZE
            MarkerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MarkerVecHC = 3 * (1 - ((MarkerVecHC - minThick_HC) / maxThick_HC))

            MarkerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MarkerVecTAU = 3 * (1 - ((MarkerVecTAU - minThick_TAU) / maxThick_TAU))

            MarkerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MarkerVecTDP = 3 * (1 - ((MarkerVecTDP - minThick_TDP) / maxThick_TDP))

            # Clear figure
            plt.clf()
            
            thicknessAtPath3DMapping(covMatList_W_Vol, NetworkDataGeneral, CoM_HC, CoM_TAU, CoM_TDP, LabelNames_Path_HC, LabelNames_Path_TAU, LabelNames_Path_TDP, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVecHC, colorVecTAU, colorVecTDP, outputDir, pthresh, fd_val, FD = True, W_Score = True)
