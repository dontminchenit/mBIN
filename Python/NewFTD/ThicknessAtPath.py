import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys
from scipy.stats import norm
from PIL import Image
from sklearn.linear_model import LinearRegression

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from analysisVisualization import drawScatterplot, nonZeroDegCmp, nonZeroDegCorr, drawCovMatrix, cropImg4, cropImg2, cropImg6
from plotNetwork3 import plotNetwork3 
from plotNetwork3_Individual import plotNetwork3Individual 
from FTDHelperFunctions import thicknessCovMatGenerate, sexToBinary, generateZScore,generateWScore, drawthicknessboxplot, save3DGraph, fixedDensity

def thicknessAtPath(outputDir, NetworkDataGeneral, pathCoM, LabelNames_Path, TAUthicknessAtPath, TDPthicknessAtPath, HCthicknessAtPath, TAUvolumeAtPath, TDPvolumeAtPath, HCvolumeAtPath, HCResults, PatientTAUResults, PatientTDPResults, TAU_missing_index_GM, TDP_missing_index_GM, sn, plotON = True):
    # Covert (20, 3, 2) --> (40, 3) / first 20 rows {L} and last 20 rows {R} / Order same as pathNames_3D_Map repeated 2 times
    pathCoM = np.vstack((pathCoM[:, :, 0], pathCoM[:, :, 1]))

    # Convert TAUthicknessAtPath, TDPthicknessAtPath, HCthicknessAtPath to appropriate shape
    # 26 x 20 x 2 --> 26 x 40 / 30 x 20 x 2 --> 30 x 40 / 54 x 20 x 2 --> 54 x 40 / Where 40 is IN ORDER TO pathT Columns {L, R}
    HCthicknessAtPath = np.hstack((HCthicknessAtPath[:, :, 0], HCthicknessAtPath[:, :, 1]))# 54 x 40
    TAUthicknessAtPath = np.hstack((TAUthicknessAtPath[:, :, 0], TAUthicknessAtPath[:, :, 1])) # 26 x 40
    TDPthicknessAtPath = np.hstack((TDPthicknessAtPath[:, :, 0], TDPthicknessAtPath[:, :, 1])) # 30 x 40

    HCvolumeAtPath = np.hstack((HCvolumeAtPath[:, :, 0], HCvolumeAtPath[:, :, 1]))# 54 x 40
    TAUvolumeAtPath = np.hstack((TAUvolumeAtPath[:, :, 0], TAUvolumeAtPath[:, :, 1])) # 26 x 40
    TDPvolumeAtPath = np.hstack((TDPvolumeAtPath[:, :, 0], TDPvolumeAtPath[:, :, 1])) # 30 x 40

    # # Drop the Regions where the Pathology Data is Missing !!
    # TAUthicknessAtPath = np.delete(TAUthicknessAtPath, TAU_missing_index_GM, axis = 1)
    # TDPthicknessAtPath = np.delete(TDPthicknessAtPath, TDP_missing_index_GM, axis = 1)

    # TAUvolumeAtPath = np.delete(TAUvolumeAtPath, TAU_missing_index_GM, axis = 1)
    # TDPvolumeAtPath = np.delete(TDPvolumeAtPath, TDP_missing_index_GM, axis = 1)    
    
    # Sanity Check
    assert(HCthicknessAtPath.shape == (54, 40))
    assert(TAUthicknessAtPath.shape == (26, 40))
    assert(TDPthicknessAtPath.shape == (30, 40))
    assert(HCvolumeAtPath.shape == (54, 40))
    assert(TAUvolumeAtPath.shape == (26, 40))
    assert(TDPvolumeAtPath.shape == (30, 40))

    # N = 40 (20 x 2)
    N = sn * 2
    
    ################### Z SCORE ###################
    # Get Mean/STD of the thickness values for HC / Mean, STD of HC for each 40 regions / shape: (40,)
    hc_Mean_Thick = np.nanmean(HCthicknessAtPath, axis=0)
    hc_SD_Thick = np.nanstd(HCthicknessAtPath, axis=0, ddof=0) # ddof parameter is set to 0, which specifies that the divisor should be N, where N is the number of non-NaN elements in the array
    
    # Get Z score for thickness values
    HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z = generateZScore(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, hc_Mean_Thick, hc_SD_Thick)

    # Get Mean/STD of the VOLUME values for HC / Mean, STD of HC for each 400 regions / 1 x 400
    hc_Mean_Vol = np.nanmean(HCvolumeAtPath, axis=0)
    hc_SD_Vol = np.nanstd(HCvolumeAtPath, axis=0, ddof=0) 

    # Get Z score for VOLUME values
    HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol = generateZScore(HCvolumeAtPath, TAUvolumeAtPath, TDPvolumeAtPath, hc_Mean_Vol, hc_SD_Vol)
    
    # Sanity Check
    assert(HCthicknessAtPath_z.shape == HCthicknessAtPath.shape)
    assert(TAUthicknessAtPath_z.shape == TAUthicknessAtPath.shape)
    assert(TDPthicknessAtPath_z.shape == TDPthicknessAtPath.shape)
    assert(HCthicknessAtPath_z_Vol.shape == HCthicknessAtPath.shape)
    assert(TAUthicknessAtPath_z_Vol.shape == TAUthicknessAtPath.shape)
    assert(TDPthicknessAtPath_z_Vol.shape == TDPthicknessAtPath.shape)

    ################### W SCORE ###################
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
    HC_w_Vol, TAU_w_Vol, TDP_w_Vol = generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCvolumeAtPath, TAUvolumeAtPath, TDPvolumeAtPath)

    # Sanity Check
    assert(HC_w.shape == HCthicknessAtPath.shape)
    assert(TAU_w.shape == TAUthicknessAtPath.shape)
    assert(TDP_w.shape == TDPthicknessAtPath.shape)
    assert(HC_w_Vol.shape == HCthicknessAtPath.shape)
    assert(TAU_w_Vol.shape == TAUthicknessAtPath.shape)
    assert(TDP_w_Vol.shape == TDPthicknessAtPath.shape)

    # Plotting and Saving figure
    saveName = 'At_Path'

    # Generate Box plot of Thickness/Z-Score/W-Score Distributions (Mean across subjects)
    drawthicknessboxplot(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, ['HC', 'TAU', 'TDP'], 'Mean Thickness', saveName, outputDir, f"FTD_Distribution_{saveName}_Original.tif")
    drawthicknessboxplot(HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, ['HC', 'TAU', 'TDP'], 'Mean Z-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_ZScore.tif")
    drawthicknessboxplot(HC_w, TAU_w, TDP_w, ['HC', 'TAU', 'TDP'], 'Mean W-Score', saveName, outputDir, f"FTD_Distribution_{saveName}_WScore.tif")

    # P value threshold
    # pthresh_list = [0.05, 0.01, 0.005, 0.001]
    pthresh_list = [0.05]
    for x in range(len(pthresh_list)):
         # p value threshold
        pthresh = pthresh_list[x]
        cov_thresh = 0.1 # originally 0.1 --> But this makes covMatHC all NaN
        
        # COMPUTING CovMat [THICKNESS]
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU
        covMatList_Original = thicknessCovMatGenerate(N, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Z = thicknessCovMatGenerate(N, HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / W score VERSION!!
        covMatList_W = thicknessCovMatGenerate(N, HC_w, TAU_w, TDP_w, pthresh, cov_thresh)

        # COMPUTING CovMat [VOLUME]
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Original_Vol = thicknessCovMatGenerate(N, HCvolumeAtPath, TAUvolumeAtPath, TDPvolumeAtPath, pthresh, cov_thresh)
        # covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU / Z score VERSION!!
        covMatList_Z_Vol = thicknessCovMatGenerate(N, HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol, pthresh, cov_thresh)
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

        
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / ORIGINAL VALUE
        ################ THICKNESS ################
        nonZeroDegCorr(HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, covMatList_Original[0], covMatList_Original[1], covMatList_Original[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original.tif", linear_regression = True)
        # Clear figure
        plt.clf()
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / Z SCORE
        nonZeroDegCorr(HCthicknessAtPath_z, TAUthicknessAtPath_z, TDPthicknessAtPath_z, covMatList_Z[0], covMatList_Z[1], covMatList_Z[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score.tif", linear_regression = True)
        # Clear figure
        plt.clf()
        # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd / Z SCORE
        nonZeroDegCorr(HC_w, TAU_w, TDP_w, covMatList_W[0], covMatList_W[1], covMatList_W[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Thickness', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score.tif", linear_regression = True)
        # Clear figure
        plt.clf()

        ################ VOLUME ################
        nonZeroDegCorr(HCvolumeAtPath, TAUvolumeAtPath, TDPvolumeAtPath, covMatList_Original_Vol[0], covMatList_Original_Vol[1], covMatList_Original_Vol[2], 'FTD_HC_' + saveName, 'FTD_TAU_' + saveName, 'FTD_TDP_' + saveName, 'Degree', 'Mean Volume', outputDir, f"FTD_nonZero_degCorr_{saveName}_Original_VOLUME.tif", linear_regression = True)
        plt.clf()
        nonZeroDegCorr(HCthicknessAtPath_z_Vol, TAUthicknessAtPath_z_Vol, TDPthicknessAtPath_z_Vol, covMatList_Z_Vol[0], covMatList_Z_Vol[1], covMatList_Z_Vol[2], 'FTD_HC_' + saveName + '_Z_Score', 'FTD_TAU_' + saveName + '_Z_Score', 'FTD_TDP_' + saveName + '_Z_Score', 'Degree', 'Mean Z-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_Z_Score_VOLUME.tif", linear_regression = True)
        plt.clf()
        nonZeroDegCorr(HC_w_Vol, TAU_w_Vol, TDP_w_Vol, covMatList_W_Vol[0], covMatList_W_Vol[1], covMatList_W_Vol[2], 'FTD_HC_' + saveName + '_W_Score', 'FTD_TAU_' + saveName + '_W_Score', 'FTD_TDP_' + saveName + '_W_Score', 'Degree', 'Mean W-Score (Volume)', outputDir, f"FTD_nonZero_degCorr_{saveName}_W_Score_VOLUME.tif", linear_regression = True)
        plt.clf()

        # 5) Plot various Covaraiance Matrix 
        if plotON:

            plotCovMat = False

            if plotCovMat: # Draw various Covaraiance Matrix
                # Covariance Matrix of Control, TAU, TDP
                drawCovMatrix(covMatList_Original[0], LabelNames_Path, LabelNames_Path, 'FTD HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_HC.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[1], LabelNames_Path, LabelNames_Path, 'FTD TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[2], LabelNames_Path, LabelNames_Path, 'FTD TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP.png", annot_bool = False)

                # Covariance Matrix - TAU, TDP comparing to Control
                drawCovMatrix(covMatList_Original[3], LabelNames_Path, LabelNames_Path, 'FTD TAU > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_HC_pthresh_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[4], LabelNames_Path, LabelNames_Path, 'FTD TDP > HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_HC_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[5], LabelNames_Path, LabelNames_Path, 'FTD TAU < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_lt_HC_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[6], LabelNames_Path, LabelNames_Path, 'FTD TDP < HC Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_lt_HC_{pthresh}.png", annot_bool = False)

                # Covariance Matrix  - comparing TAU vs TDP
                drawCovMatrix(covMatList_Original[7], LabelNames_Path, LabelNames_Path, 'FTD TAU > TDP Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TAU_gt_TDP_{pthresh}.png", annot_bool = False)
                drawCovMatrix(covMatList_Original[8], LabelNames_Path, LabelNames_Path, 'FTD TDP > TAU Covariance Matrix', outputDir, f"FTD_Thickness_{saveName}_TDP_gt_TAU_{pthresh}.png", annot_bool = False)
            

            # 6) Plot Covariance Matrix (Network) onto 3D Atlas / For Original Values (because Z-Score generates the same figure)
            cRange = [0, 1]
            
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

            # thick_HC --> N x 40
            # np.nanmean(thick_HC_exp, axis=0) --> 1 x 40
            
            # MakerVec for HC, TAU and TDP -> FOR NODE SIZE
            MakerVecHC = np.nanmean(thick_HC_exp, axis=0)
            MakerVecHC = 3 * (1 - ((MakerVecHC - minThick_HC) / maxThick_HC))

            MakerVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            MakerVecTAU = 3 * (1 - ((MakerVecTAU - minThick_TAU) / maxThick_TAU))

            MakerVecTDP = np.nanmean(thick_TDP_exp, axis=0)
            MakerVecTDP = 3 * (1 - ((MakerVecTDP - minThick_TDP) / maxThick_TDP))

            # # colorVec for YesAD and NoAD
            # colorVecTAU = np.nanmean(thick_TAU_exp, axis=0)
            # colorVecTDP = np.nanmean(thick_TDP_exp, axis=0)

            # Node color --> Set as red (because cm.jet: 1 --> Red)
            # colorVecTAU = np.ones(N)
            # colorVecTDP = np.ones(N)
            colorVec = np.ones(N)

            # Originally 0 - but using 3 bc 0 is not yet implemented
            displayType = 3

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

            covMatNamelist = ['covMatHC_Path', 'covMatTAU_Path', 'covMatTDP_Path', 'cmpCovTAU_gt_HC_Path', 'cmpCovTDP_gt_HC_Path', 'cmpCovTAU_lt_HC_Path', 'cmpCovTDP_lt_HC_Path', 'cmpCovTAU_gt_TDP_Path', 'cmpCovTDP_gt_TAU_Path']
            covMatNamelist = ['covMatHC_Path', 'covMatTAU_Path', 'covMatTDP_Path', 'cmpCovTAU_gt_HC_Path', 'cmpCovTDP_gt_HC_Path', 'cmpCovTAU_lt_HC_Path', 'cmpCovTDP_lt_HC_Path', 'cmpCovTAU_gt_TDP_Path', 'cmpCovTDP_gt_TAU_Path', f'cmpCovTAU_gt_HC_FD_{fd_val}_Path', f'cmpCovTDP_gt_HC_FD_{fd_val}_Path', f'cmpCovTAU_lt_HC_FD_{fd_val}_Path', f'cmpCovTDP_lt_HC_FD_{fd_val}_Path', f'cmpCovTAU_gt_TDP_FD_{fd_val}_Path', f'cmpCovTDP_gt_TAU_FD_{fd_val}_Path']
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

                #### Modify CovMat - To match the Missing Pathology Dataset (NaN Values) ####
                if j == 0:
                    currCovMat = covMatlist[j]
                    currPathCoM = pathCoM
                    currLabelNames_Path = LabelNames_Path
                    currMakerVec = MakerVec
                    currColorVec = colorVec
                elif j == 1 or j == 3 or j == 5 or j == 7: # TAU part
                    currCovMat = np.delete(np.delete(covMatlist[j], TAU_missing_index_GM, axis = 0), TAU_missing_index_GM, axis = 1) # (40 x 40) --> (35 x 35)
                    currPathCoM = np.delete(pathCoM, TAU_missing_index_GM, axis = 0) #(40, 3) --> (35, 3)
                    currLabelNames_Path = np.delete(LabelNames_Path, TAU_missing_index_GM, axis=0) # (40,) --> (35,)
                    currMakerVec = np.delete(MakerVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                    currColorVec = np.delete(colorVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                else: # TDP part
                    currCovMat = np.delete(np.delete(covMatlist[j], TDP_missing_index_GM, axis = 0), TDP_missing_index_GM, axis = 1) # (40 x 40) --> (34 x 34)
                    currPathCoM = np.delete(pathCoM, TDP_missing_index_GM, axis = 0) #(40, 3) --> (34, 3)
                    currLabelNames_Path = np.delete(LabelNames_Path, TDP_missing_index_GM, axis=0) # (40,) --> (34,)
                    currMakerVec = np.delete(MakerVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)
                    currColorVec = np.delete(colorVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)

                # [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]
                plotNetwork3Individual(fig_atlas, currCovMat, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currPathCoM, currLabelNames_Path, cRange, currMakerVec, currColorVec, 3, showLabels = 0, covType = covType)

                fig_atlas.tight_layout()

                # Save Graph Object
                save3DGraph(currCovMat, outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}')

                # Save the figure
                plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)', dpi=1000, bbox_inches='tight')

                # Read figure to crop
                imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)' + '.png')

            thickness_comb_img = cropImg6(imglist[3:9])
            thickness_comb_FD_img = cropImg6(imglist[9:])

            thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(Original).png')
            thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(FixedDensity).png')

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

            # FOR NODE SIZE (Normalized Volume)
            thick_HC_exp = HCvolumeAtPath
            thick_TAU_exp = TAUvolumeAtPath
            thick_TDP_exp = TDPvolumeAtPath

            # Get the MAX/MIN Thickness values
            minThick_HC = np.nanmin(np.nanmean(thick_HC_exp, axis=0))
            minThick_TAU = np.nanmin(np.nanmean(thick_TAU_exp, axis=0))
            minThick_TDP = np.nanmin(np.nanmean(thick_TDP_exp, axis=0))

            vanishing_val = 0.2
            maxThick_HC = np.nanmax(np.nanmean(thick_HC_exp, axis=0) - minThick_HC) + vanishing_val
            maxThick_TAU = np.nanmax(np.nanmean(thick_TAU_exp, axis=0) - minThick_TAU) + vanishing_val
            maxThick_TDP = np.nanmax(np.nanmean(thick_TDP_exp, axis=0) - minThick_TDP) + vanishing_val

            # thick_HC --> N x 40
            # np.nanmean(thick_HC_exp, axis=0) --> 1 x 40
            
            # MakerVec for HC, TAU and TDP -> FOR NODE SIZE
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

                #### Modify CovMat - To match the Missing Pathology Dataset (NaN Values) ####
                if j == 0:
                    currCovMat = covMatlist[j]
                    currPathCoM = pathCoM
                    currLabelNames_Path = LabelNames_Path
                    currMakerVec = MakerVec
                    currColorVec = colorVec
                elif j == 1 or j == 3 or j == 5 or j == 7: # TAU part
                    currCovMat = np.delete(np.delete(covMatlist[j], TAU_missing_index_GM, axis = 0), TAU_missing_index_GM, axis = 1) # (40 x 40) --> (35 x 35)
                    currPathCoM = np.delete(pathCoM, TAU_missing_index_GM, axis = 0) #(40, 3) --> (35, 3)
                    currLabelNames_Path = np.delete(LabelNames_Path, TAU_missing_index_GM, axis=0) # (40,) --> (35,)
                    currMakerVec = np.delete(MakerVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                    currColorVec = np.delete(colorVec, TAU_missing_index_GM, axis = 0) # (40,) --> (35,)
                else: # TDP part
                    currCovMat = np.delete(np.delete(covMatlist[j], TDP_missing_index_GM, axis = 0), TDP_missing_index_GM, axis = 1) # (40 x 40) --> (34 x 34)
                    currPathCoM = np.delete(pathCoM, TDP_missing_index_GM, axis = 0) #(40, 3) --> (34, 3)
                    currLabelNames_Path = np.delete(LabelNames_Path, TDP_missing_index_GM, axis=0) # (40,) --> (34,)
                    currMakerVec = np.delete(MakerVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)
                    currColorVec = np.delete(colorVec, TDP_missing_index_GM, axis = 0) # (40,) --> (34,)

                # [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]
                plotNetwork3Individual(fig_atlas, currCovMat, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currPathCoM, currLabelNames_Path, cRange, currMakerVec, currColorVec, 3, showLabels = 0, covType = covType)

                fig_atlas.tight_layout()

                # Save Graph Object
                save3DGraph(currCovMat, outputDir, covMatNamelist[j] + f'_Path_pthresh_{str(pthresh).split(".")[1]}_WScore')

                # Save the figure
                plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)', dpi=1000, bbox_inches='tight')

                # Read figure to crop
                imglist_W[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)' + '.png')

            thickness_comb_img_W = cropImg6(imglist_W[3:9])
            thickness_comb_FD_img_W = cropImg6(imglist_W[9:])

            thickness_comb_img_W.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(WSCORE).png')
            thickness_comb_FD_img_W.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(WSCORE_FixedDensity).png')
