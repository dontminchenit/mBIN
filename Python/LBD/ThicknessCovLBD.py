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

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from analysisVisualization import drawScatterplot, nonZeroDegCmp, nonZeroDegCorr, drawCovMatrix
from plotNetwork3 import plotNetwork3 

# Helper Functions
def test_corr_sig(corr1, corr2, ncorr1, ncorr2):
    z_corr1 = np.arctanh(corr1)
    z_corr2 = np.arctanh(corr2)
    z_obs = (z_corr1 - z_corr2) / ((1 / (ncorr1 - 3)) + (1 / (ncorr2 - 3))) ** 0.5
    probs = norm.cdf(z_obs)
    return probs, z_obs

import numpy as np

def efficiency_bin(A, local=False):
    """
    Calculate efficiency based on input data A.
    
    Parameters:
    -----------
    A : numpy.ndarray
        Input data array.
    local : bool, optional
        Flag to specify whether to calculate local efficiency (default is False).
    
    Returns:
    --------
    E : numpy.ndarray
        Efficiency values.
    """
    n = len(A)                          # number of nodes
    np.fill_diagonal(A, 0)              # clear diagonal
    A = np.double(A != 0)               # enforce double precision

    if local:                           # local efficiency
        E = np.zeros(n)
        for u in range(n):
            V = np.nonzero(A[u,:] | A[:,u])[0]       # neighbors
            sa = A[u,V] + A[V,u].T                  # symmetrized adjacency vector
            e = distance_inv(A[V,:][:,V])           # inverse distance matrix
            se = e + e.T                            # symmetrized inverse distance matrix
            numer = np.sum(np.sum(sa.T * sa) * se) / 2  # numerator
            if numer != 0:
                denom = np.sum(sa) ** 2 - np.sum(sa ** 2)  # denominator
                E[u] = numer / denom                  # local efficiency
    else:                                         # global efficiency
        e = distance_inv(A)
        E = np.sum(e) / (n ** 2 - n)
    
    return E


def distance_inv(A_):
    l = 1                                             # path length
    Lpath = A_                                        # matrix of paths l
    D = A_.copy()                                      # distance matrix
    n_ = len(A_)

    Idx = True
    while np.any(Idx):
        l += 1
        Lpath = Lpath @ A_
        Idx = (Lpath != 0) & (D == 0)
        D[Idx] = l

    D[~D | np.eye(n_).astype(bool)] = np.inf         # assign inf to disconnected nodes and to diagonal
    D = 1. / D                                       # invert distance

    return D


def ThicknessCovLBD(NetworkDataGeneral, pathLUT, CtrlResults, AllResults, demoT, outputDir, measuresT_matched, plotON = True):
    # List of network names
    networkNames=['DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All']

    # Define figure - subplot with 2 rows and 4 columns
    fig = plt.figure()

    # 1000 rows for now
    NN = 1000

    # Create a Pandas Dataframe of 1000 rows x 13 columns (Initally all empty - NaN)
    NetworkNames = [np.nan for _ in range(NN)]
    ROIName1 = [np.nan for _ in range(NN)]
    ROIName1_Path = [np.nan for _ in range(NN)]
    ROIName2 = [np.nan for _ in range(NN)]
    ROIName2_Path = [np.nan for _ in range(NN)]
    r_YesAD = [np.nan for _ in range(NN)]
    r_NoAD = [np.nan for _ in range(NN)]
    N_YesAD = [np.nan for _ in range(NN)]
    N_NoAD = [np.nan for _ in range(NN)]
    z_YesAD_gt_NoAD = [np.nan for _ in range(NN)]
    z_NoAD_gt_YesAD = [np.nan for _ in range(NN)]
    p_YesAD_gt_NoAD = [np.nan for _ in range(NN)]
    p_NoAD_gt_YesAD = [np.nan for _ in range(NN)]

    SigEdgesT = pd.DataFrame({'NetworkNames': NetworkNames, 'ROIName1': ROIName1, 'ROIName1_Path': ROIName1_Path, 'ROIName2': ROIName2, 'ROIName2_Path': ROIName2_Path, 'r_YesAD': r_YesAD, 'r_NoAD': r_NoAD, 'N_YesAD': N_YesAD, 'N_NoAD': N_NoAD, 'z_YesAD_gt_NoAD': z_YesAD_gt_NoAD, 'z_NoAD_gt_YesAD': z_NoAD_gt_YesAD, 'p_YesAD_gt_NoAD': p_YesAD_gt_NoAD, 'p_NoAD_gt_YesAD': p_NoAD_gt_YesAD})
    SigEdgesCounter = 0

    # Denothing what NetworkNames to use
    # 5 --> 'Default'
    ii = 5
    
    ###################################################################################################################################################
    # Label_Name2_list: Names of 400 Regions in the Schaffer Atlas (in order with numlab = 400)
    # currNetwork: List of Boolean denoting if Label_Name2_list contains networkNames('Default') or not / Boolean list of length 400
    # pathLUTSorted: Sort the pathology LUT by Label_ID400x7 1 to 400 / LUT matching 400 thickness regions to PathName
    # pathNames: PathNames that contains 'Default' as list (in order of the 400 regions)
    # regNamesRaw: Label_Name2 that contains current Network ('Default') as list (in order of the 400 regions)
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

    # saveName --> this specific case: 'Default'
    saveName = networkNames[ii - 1] # -1 added because of index difference between python and matlab

    # Sort the pathology LUT by Label_ID400x7
    # pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))
    pathLUTSorted = pathLUT.sort_values(by=['Label_ID400x7']) # Sort by label 1 to 400
    
    # Get the PathNames of the current Network ('Default') inside the pathLUTSorted dataset as list / Out of 400 Regions there are 91 regions that contains 'Default'
    pathNames = pathLUTSorted[currNetwork]['PathName'].tolist()

    # Get the Label_Name2 that contains current Network ('Default') as list
    regNamesRaw = list(compress(Label_Name2_list, currNetwork))

    # Length of regNamesRaw (=91)
    SN = len(regNamesRaw)

    # Change the regNamesRaw (_ to -)
    regNames = []

    for j in range(SN):
        nameSplit = regNamesRaw[j].replace('_', '-')
        regNames.append(nameSplit)

    # Get current Center of Mass for current Network / 91
    currCoM = (NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['CoM'][0, 0])[currNetwork]

    # N = sum(currNetwork) / basically same as SN / 91
    N = SN

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CtrlResults: thickness mean and volume total values for [Control MRI data IDs x lables (numLab)]
    # AllResults: thickness mean and volume total values for [Patient MRI data IDs x lables (numLab = 400 regions in the sch region)]
    # demoT: MRI Data - Demographics 

    # Get thickness mean values from Control MRI Data for only current Network [Control MRI data IDs x N (91)] --> 169 Control x 91 regions
    thickCtrl = CtrlResults['Thickness']['Mean'][:,currNetwork]

    # Get thickness mean values from Patient MRI Data for only current Network [Patient MRI data IDs (with or without AD) x N (91)] --> 75 patients x 91 regions
    # But further divide them into with AD or without AD
    thickyesAD = AllResults['Thickness']['Mean'][:, currNetwork] # With AD / 33 x 91
    thickyesAD = thickyesAD[demoT['LBD0_LBDAD1'] == 1, :]

    thicknoAD = AllResults['Thickness']['Mean'][:, currNetwork] # Without AD / 42 x 91
    thicknoAD = thicknoAD[demoT['LBD0_LBDAD1'] == 0, :]

    # Get Mean/STD of the thickness values for Control / Mean, STD of Controls for each 91 regions / 1 x 91
    Ctrl_Mean = np.nanmean(thickCtrl, axis=0)
    Ctrl_SD = np.nanstd(thickCtrl, axis=0, ddof=0) # ddof parameter is set to 0, which specifies that the divisor should be N, where N is the number of non-NaN elements in the array

    # YesAD_W: W score of Thickness values for Patient Dataset with AD
    YesAD_W = np.empty(thickyesAD.shape)
    for i in range(thickyesAD.shape[0]):
        YesAD_W[i, :] = (thickyesAD[i, :] - Ctrl_Mean) / Ctrl_SD

    # NoAD_W: W score of Thickness values for Patient Dataset without AD
    NoAD_W = np.empty(thicknoAD.shape)
    for i in range(thicknoAD.shape[0]):
        NoAD_W[i, :] = (thicknoAD[i, :] - Ctrl_Mean) / Ctrl_SD

    # Ctrl_W: W score of Thickness values for Control Dataset [QUESTION - [Why set to thicknoAD shape??]]
    Ctrl_W = np.empty(thicknoAD.shape)
    for i in range(thicknoAD.shape[0]):
        Ctrl_W[i, :] = (thickCtrl[i, :] - Ctrl_Mean) / Ctrl_SD

    # # Mean of W scores for Patient Dataset with or without AD
    mean_YesAD_w = np.nanmean(YesAD_W, axis=0) # 1 x 91 / 33 Patients
    mean_NoAD_w = np.nanmean(NoAD_W, axis=0) # 1 x 91 / 42 Patients

    # Covariance Matrix - Initialization (N = number of regions in currNetwork = 91)
    covMatCtrl = np.full((N,N), np.nan)
    covMatyesAD = np.full((N,N), np.nan)
    covMatnoAD = np.full((N,N), np.nan)

    cmpCovYes_gt_Ctrl = np.full((N,N), np.nan)
    cmpCovNo_gt_Ctrl = np.full((N,N), np.nan)

    cmpCovYes_lt_Ctrl = np.full((N,N), np.nan)
    cmpCovNo_lt_Ctrl = np.full((N,N), np.nan)

    cmpCovYes_gt_No = np.full((N,N), np.nan)
    cmpCovNo_gt_Yes = np.full((N,N), np.nan)

    # Calculating Covariance Matrix
    # P value threshold
    pthresh = 0.01

    for i in range(N):
        for j in range(N):

            if i != j: # When it's not the same region in the current Network
                NCtrl = np.sum(~np.isnan(thickCtrl[:, i]) & ~np.isnan(thickCtrl[:, j])) # Get the sum of regions that are not NaN
                NYes = np.sum(~np.isnan(thickyesAD[:, i]) & ~np.isnan(thickyesAD[:, j])) # Get the sum of regions that are not NaN
                NNo = np.sum(~np.isnan(thicknoAD[:, i]) & ~np.isnan(thicknoAD[:, j])) # Get the sum of regions that are not NaN

                if NCtrl > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatCtrl[i, j] = ma.corrcoef(ma.masked_invalid(thickCtrl[:,i]), ma.masked_invalid(thickCtrl[:,j]))[0, 1]
                    if covMatCtrl[i, j] < 0.1:
                        covMatCtrl[i, j] = np.nan
                        NCtrl = 0 # reset for below case

                if NYes > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatyesAD[i, j] = ma.corrcoef(ma.masked_invalid(thickyesAD[:,i]), ma.masked_invalid(thickyesAD[:,j]))[0, 1]
                    if covMatyesAD[i, j] < 0.1:
                        covMatyesAD[i, j] = np.nan

                if NNo > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatnoAD[i, j] = ma.corrcoef(ma.masked_invalid(thicknoAD[:,i]), ma.masked_invalid(thicknoAD[:,j]))[0, 1]
                    if covMatnoAD[i, j] < 0.1:
                        covMatnoAD[i, j] = np.nan

                if NCtrl > 3 and NYes > 3: # NOT REACHABLE???
                    cmpCovYes_gt_Ctrl[i, j] = test_corr_sig(covMatCtrl[i, j], covMatyesAD[i, j], NCtrl, NYes)[0] < pthresh # get only probs
                    cmpCovYes_lt_Ctrl[i, j] = test_corr_sig(covMatyesAD[i, j], covMatCtrl[i, j], NYes, NCtrl)[0] < pthresh # get only probs

                if NCtrl > 3 and NNo > 3: # NOT REACHABLE???
                    cmpCovNo_gt_Ctrl[i, j] = test_corr_sig(covMatCtrl[i, j], covMatnoAD[i, j], NCtrl, NNo)[0] < pthresh # get only probs
                    cmpCovNo_lt_Ctrl[i, j] = test_corr_sig(covMatnoAD[i, j], covMatCtrl[i, j], NNo, NCtrl)[0] < pthresh # get only probs

                if NYes > 3 and NNo > 3:
                    cmpCovYes_gt_No[i, j] = test_corr_sig(covMatnoAD[i, j], covMatyesAD[i, j], NNo, NYes)[0] < pthresh # get only probs
                    cmpCovNo_gt_Yes[i, j] = test_corr_sig(covMatyesAD[i, j], covMatnoAD[i, j], NYes, NNo)[0] < pthresh # get only probs
                
                # Fill in the SigEdgesT pandas dataframe we created earlier
                if (cmpCovYes_gt_No[i, j] == 1 or cmpCovNo_gt_Yes[i, j] == 1) and i < j: # record significant edges
                    SigEdgesCounter += 1 # increase SigEdgesCounter
                    ss = SigEdgesCounter - 1 # index difference between matlab and python

                    SigEdgesT.at[ss, 'NetworkNames'] = networkNames[ii - 1] # index difference between matlab and python
                    SigEdgesT.at[ss, 'ROIName1'] = regNames[i]
                    SigEdgesT.at[ss, 'ROIName1_Path'] = pathNames[i] # NaN is obmitted         
                    SigEdgesT.at[ss, 'ROIName2'] = regNames[j]
                    SigEdgesT.at[ss, 'ROIName2_Path'] = pathNames[j] # NaN is obmitted 
                    SigEdgesT.at[ss, 'r_YesAD'] = ma.corrcoef(ma.masked_invalid(thickyesAD[:,i]), ma.masked_invalid(thickyesAD[:,j]))[0, 1]
                    SigEdgesT.at[ss, 'r_NoAD'] = ma.corrcoef(ma.masked_invalid(thicknoAD[:,i]), ma.masked_invalid(thicknoAD[:,j]))[0, 1]
                    SigEdgesT.at[ss, 'N_YesAD'] = NYes
                    SigEdgesT.at[ss, 'N_NoAD'] = NNo

                    p, z = test_corr_sig(covMatnoAD[i, j], covMatyesAD[i, j], NNo, NYes)
                    SigEdgesT.at[ss, 'p_YesAD_gt_NoAD'] = p
                    SigEdgesT.at[ss, 'z_YesAD_gt_NoAD'] = z

                    p, z = test_corr_sig(covMatyesAD[i, j], covMatnoAD[i, j], NYes, NNo)
                    SigEdgesT.at[ss, 'p_NoAD_gt_YesAD'] = p
                    SigEdgesT.at[ss, 'z_NoAD_gt_YesAD'] = z
    
    # Save SigEdgesT as csv
    SigEdgesT.to_csv(outputDir + "/SigEdges_20230106.csv")  

    # One of the measurement you can calculate on Network (Global measurement)
    # efficiency_bin(cmpCovYes_gt_Ctrl)
    # efficiency_bin(cmpCovNo_gt_Ctrl)
    # efficiency_bin(covMatCtrl)

    # ROIyesAD ROInoAD (ROI --> restricting to the area where there is a significance diff between Control vs Patient) / Boolean 0 or 1
    ROIyes = np.sum(cmpCovYes_gt_Ctrl, axis=1, keepdims=True, where=~np.isnan(cmpCovYes_gt_Ctrl)) > 0 # 1 x 91 (regions)
    ROIno = np.sum(cmpCovNo_gt_Ctrl, axis=1, keepdims=True, where=~np.isnan(cmpCovNo_gt_Ctrl)) > 0 # 1 x 91 (regions)

    # Thickness values Mean of regions where ROIyes or ROINo. For Patient MRI Data for only current Network
    key_thickness_yes = np.mean(thickyesAD[:, ROIyes.flatten()], axis=1) # 33 x 1
    key_thickness_no = np.mean(thicknoAD[:, ROIno.flatten()], axis=1) # 42 x 1
    
    # measureT: clinical test, measurement (Traditional  Neurological tests in measuring disease severity)
    # measuresT_matched: measuresT Dataframe with only the rows with match (demoT.INDDID are in measuresT.INDDID)
    measureData = measuresT_matched['BNT  ANN CHANGE']

    measure_yes = measureData[(demoT['LBD0_LBDAD1'] == 1).tolist()].to_numpy() # with AD / 33 x 1
    measure_no = measureData[(demoT['LBD0_LBDAD1'] == 0).tolist()].to_numpy() # without AD / 42 x 1

    # Plotting and Saving figure
    
    # 1) Scatter plot for key_thickness_yes vs measure_yes
    # Ignore NaN Values
    nas = np.logical_or(np.isnan(key_thickness_yes), np.isnan(measure_yes))

    drawScatterplot(key_thickness_yes[~nas], measure_yes[~nas], 'Mean Thickness at Significant Covariance ROIs', 'MMSE Annualized Change', 
                    f"LBD+AD, {saveName}", outputDir, f"measureCmpMMSE_YES_{saveName}_.tif", linear_regression = True)
    
    # Clear figure
    plt.clf()

    # 2) Scatter plot for key_thickness_no vs measure_no
    # Ignore NaN Values
    nas = np.logical_or(np.isnan(key_thickness_no), np.isnan(measure_no))

    drawScatterplot(key_thickness_no[~nas], measure_no[~nas], 'Mean Thickness at Significant Covariance ROIs', 'MMSE Annualized Change', 
                    f"LBD-AD, {saveName}", outputDir, f"measureCmpMMSE_NO_{saveName}_.tif", linear_regression = True)
    
    # Clear figure
    plt.clf()
    
    # 3) Comparing the degree of the network between Control vs YesAD(Patient) vs NoAd(Patient)
    nonZeroDegCmp(covMatCtrl, covMatyesAD, covMatnoAD, ['Ctrl', 'Yes AD', 'No AD'], 'Degree', saveName, outputDir, f"nonZero_degCmp_{saveName}.tif")

    # Clear figure
    plt.clf()

    # 4) Plot correlation scatterplot between Degree value vs Thickness / for both YesAD and NoAd
    nonZeroDegCorr(thickyesAD, thicknoAD, covMatyesAD, covMatnoAD, saveName, 'Degree', 'Thickness', outputDir, f"nonZero_degCorr_{saveName}.tif", linear_regression = True)

    # Clear figure
    plt.clf()

    # 5) Plot various Covaraiance Matrix 
    if plotON:

        plotCovMat = False

        if plotCovMat: # Draw various Covaraiance Matrix
            # Covariance Matrix of Control, YesAD, NoAD
            drawCovMatrix(covMatCtrl, regNamesRaw, regNamesRaw, 'Controls Covariance Matrix', outputDir, f"{saveName}_ctrl.png")
            drawCovMatrix(covMatyesAD, regNamesRaw, regNamesRaw, 'LBD-yesAD Covariance Matrix', outputDir, f"{saveName}_yesAD.png")
            drawCovMatrix(covMatnoAD, regNamesRaw, regNamesRaw, 'LBD-noAD Covariance Matrix', outputDir, f"{saveName}_noAD.png")

            # Covariance Matrix - YesAD, NoAD comparing to Control
            drawCovMatrix(cmpCovYes_gt_Ctrl, regNamesRaw, regNamesRaw, 'YesAD > Ctrl Covariance Matrix', outputDir, f"{saveName}_YesGTCTRL.png")
            drawCovMatrix(cmpCovNo_gt_Ctrl, regNamesRaw, regNamesRaw, 'NoAD > Ctrl Covariance Matrix', outputDir, f"{saveName}_NOGTCTRL.png")

            # Covariance Matrix  - comparing YesAD vs NoAD
            drawCovMatrix(cmpCovYes_gt_No, regNamesRaw, regNamesRaw, 'YesAD > NoAD Covariance Matrix', outputDir, f"{saveName}_YesGTNO.png")
            drawCovMatrix(cmpCovNo_gt_Yes, regNamesRaw, regNamesRaw, 'NoAD > YesAD Covariance Matrix', outputDir, f"{saveName}_NOGTYES.png")
        

        # 6) Plot Covariance Matrix (Network) onto 3D Atlas
        cRange = [0, 1]
        
        MakerVec = np.ones(np.sum(currNetwork)) # vector all 1 of length 91

        # Copy W Score of Patients with or without AD
        thick_yesAD_exp = YesAD_W
        thick_noAD_exp = NoAD_W

        # Get the MAX/MIN Thickness values (from the YesAD and NoAD W Score)
        minThick = np.nanmin(np.concatenate([thick_yesAD_exp.flatten(), thick_noAD_exp.flatten()]))
        maxThick = np.nanmax((np.concatenate([thick_yesAD_exp.flatten(), thick_noAD_exp.flatten()])) - minThick + 0.0015)
        
        # MakerVec for YesAD and NoAD
        MakerVecYes = np.interp(np.nanmean(thick_yesAD_exp, axis=0), (np.nanmin(thick_yesAD_exp), np.nanmax(thick_yesAD_exp)), (0.35, 0.9))
        MakerVecYes = 1.0 / MakerVecYes

        MakerVecNo = np.interp(np.nanmean(thick_noAD_exp, axis=0), (np.nanmin(thick_noAD_exp), np.nanmax(thick_noAD_exp)), (0.35, 0.9))
        MakerVecNo = 1.0 / MakerVecNo

        # colorVec for YesAD and NoAD
        colorVecYes = np.nanmean(thick_yesAD_exp, axis=0)
        colorVecNo = np.nanmean(thick_noAD_exp, axis=0)

        # QUESTION - Why overwrite the values?
        colorVecYes = np.ones(np.sum(currNetwork))
        colorVecNo = np.ones(np.sum(currNetwork))

        MakerVecYes = np.ones(np.sum(currNetwork))
        MakerVecNo = np.ones(np.sum(currNetwork))

        # Originally 0 - but using 3 bc 0 is not yet implemented
        displayType = 3

        LabelNames = regNamesRaw

        # Clear figure
        plt.clf()

        # Define figure
        fig_atlas = plt.figure()
    
        pSelectCol = 1
        plotNetwork3(fig_atlas, cmpCovYes_gt_Ctrl, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectCol, colorVecYes, displayType)

        pSelectcol=2
        plotNetwork3(fig_atlas, cmpCovNo_gt_Ctrl, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecNo, pSelectcol, colorVecNo, displayType)

        pSelectcol=3
        plotNetwork3(fig_atlas, cmpCovYes_gt_No, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectcol, colorVecYes,displayType)

        pSelectcol=4
        plotNetwork3(fig_atlas, cmpCovNo_gt_Yes, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecNo, pSelectcol, colorVecNo, displayType)

        fig_atlas.tight_layout()

        # Save the figure
        plt.savefig(outputDir + f"/{saveName}_SigPlots.png", dpi=1000, bbox_inches='tight')