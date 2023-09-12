import numpy as np
import numpy.ma as ma
from scipy.stats import norm
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy import stats
import matplotlib.pyplot as plt
import os
import networkx as nx
import pickle
import sys
from PIL import Image
from itertools import compress

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from plotNetwork3_Individual import plotNetwork3Individual
from analysisVisualization import cropImg4, cropImg6

def test_corr_sig(corr1, corr2, ncorr1, ncorr2):
    """Statistical tool for calculating the corr significance

    Args:
        corr1 (float): covariance value of covMat1 [i, j]
        corr2 (float): covariance value of covMat2 [i, j]
        ncorr1 (int): number of non-Nan observations for covMat1 at [i, j]
        ncorr2 (int): number of non-Nan observations for covMat2 at [i, j]

    Returns:
        probs (float): pvalue
        z_obs (float): _description_
    """
    z_corr1 = np.arctanh(corr1)
    z_corr2 = np.arctanh(corr2)
    z_obs = (z_corr1 - z_corr2) / ((1 / (ncorr1 - 3)) + (1 / (ncorr2 - 3))) ** 0.5
    probs = norm.cdf(z_obs)
    return probs, z_obs

def save3DGraph(adjMtx, outputDir, filename):
    """save Graph as obj that we are mapping to 3D

    Args:
        adjMtx (ndarray): adjacency matrix we use to generate the Graph obj
        outputDir (str): Path location to save the analysis results
        filename (str): Name of the file we are saving
    """
    # Substitute NaN values to 0
    adjMtx[np.isnan(adjMtx)] = 0

    # Get upper triangle of adjMtx
    adjMtx = np.triu(adjMtx)

    # Get graph of adjMtx
    G = nx.Graph(adjMtx)

    # save graph object to file
    pickle.dump(G, open(outputDir + '/Graph_Objects/' + filename + '.pickle', 'wb'))

################################################################## PATHOLOGY PART ##################################################################
def pathObsThresh(TAUData, TDPData, obs_thresh):
    """Find the list Region index where there are less than pre-defined observation count (for both TAU and TDP)

    Args:
        TAUData (ndarray): log %AO of TAU Pathology Data
        TDPData (ndarray): log %AO of TDP Pathology Data
        obs_thresh (int): PreDefined observation count. Any regions' observation less than this value would be dicarded
    
    Returns:
        TAU_missing_index (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TAU)
        TDP_missing_index (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count (for TDP)
    """
     # Get list containing number of observed subjects per region (list of length 40) / TAU and TDP
    obs_TAU = []
    obs_TDP = []

    for t in range(TAUData.shape[1]): # iterate over the rows of the 2d array / 40
        non_nan_count = np.count_nonzero(~np.isnan(TAUData[:, t])) # Number of Non-NaN in this column
        obs_TAU.append(non_nan_count)
    
    for t in range(TDPData.shape[1]):
        non_nan_count = np.count_nonzero(~np.isnan(TDPData[:, t])) # Number of Non-NaN in this column
        obs_TDP.append(non_nan_count)

    # TO numpy array
    obs_TAU = np.array(obs_TAU)
    obs_TDP = np.array(obs_TDP)   
    
    TAU_missing_index = np.argwhere(obs_TAU < obs_thresh).flatten()
    TDP_missing_index = np.argwhere(obs_TDP < obs_thresh).flatten()

    return TAU_missing_index, TDP_missing_index

def pathCovGen(TAUData, TDPData, pthresh, cov_thresh):
    """Function that calculates covariance matrix(s) for pathology dataset.
    4 Cov Mat
    1) TAU
    2) TDP
    3) TAU_gt_TDP
    4) TDP_gt_TAU
    5) covTAU_gt_TDP_raw
    6) covTDP_gt_TAU_raw

    Args:
        TAUData (ndarray): log %AO of TAU Pathology Data
        TDPData (ndarray): log %AO of TDP Pathology Data
        pthresh (float): p-value threshold
        cov_thresh (float): covariance matrix threshold (just to make suere there is no noise)

    Returns:
        [covTAU, covTDP, covTAU_gt_TDP, covTDP_gt_TAU, covTAU_gt_TDP_raw, covTDP_gt_TAU_raw] (list of ndarray): list of 6 cov mat we calculated
    """
    # Sanity Check
    assert(TAUData.shape[1] == TDPData.shape[1])

    # Define N (Number of regions)
    N = TAUData.shape[1]

    # Covariance matrix for TAU
    covTAU= np.full((N,N), np.nan)

    # Covariance matrix for TDP
    covTDP= np.full((N,N), np.nan)

    # For Significance difference in Covariance value
    covTAU_gt_TDP = np.full((N,N), np.nan)
    covTDP_gt_TAU = np.full((N,N), np.nan)

    covTAU_gt_TDP_raw = np.full((N,N), np.nan)
    covTDP_gt_TAU_raw = np.full((N,N), np.nan)

    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of log %AO for two different anatomical regions for TAU/TDP both not NaN
                nTAU = np.sum(~np.isnan(TAUData[:,i]) & ~np.isnan(TAUData[:,j]))
                nTDP = np.sum(~np.isnan(TDPData[:,i]) & ~np.isnan(TDPData[:,j]))
                
                if nTAU > 3:
                    covTAU[i, j] = ma.corrcoef(ma.masked_invalid(TAUData[:,i]), ma.masked_invalid(TAUData[:,j]))[0, 1]
                    if covTAU[i, j] < cov_thresh:
                        covTAU[i, j] = np.nan
                        nTAU = 0

                if nTDP > 3:
                    covTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covTDP[i, j] < cov_thresh:
                        covTDP[i, j] = np.nan
                        nTDP = 0

                if nTAU > 3 and nTDP >3:
                    covTAU_gt_TDP[i,j] = test_corr_sig(covTDP[i,j],covTAU[i,j],nTDP,nTAU)[0] < pthresh 
                    covTAU_gt_TDP_raw[i,j] = test_corr_sig(covTDP[i,j],covTAU[i,j],nTDP,nTAU)[0] # Get sig value

                    covTDP_gt_TAU[i,j] = test_corr_sig(covTAU[i,j],covTDP[i,j],nTAU,nTDP)[0] < pthresh
                    covTDP_gt_TAU_raw[i,j] = test_corr_sig(covTAU[i,j],covTDP[i,j],nTAU,nTDP)[0] # get sig value

    return [covTAU, covTDP, covTAU_gt_TDP, covTDP_gt_TAU, covTAU_gt_TDP_raw, covTDP_gt_TAU_raw]

def path3DMapping(covMatlist, NetworkDataGeneral, CoM_TAU, pathNames_TAU, MarkerVecTAU, colorVecTAU, CoM_TDP, pathNames_TDP, MarkerVecTDP, colorVecTDP, cRange, outputDir, suffix_M, fd_val):
    """3D mapping of the Pathology Data Analysis

    Args:
        covMatlist (list of ndarray): list of 4 cov mat we calculated for pathology dataset
        NetworkDataGeneral (obj): Preconstructed atlas data
        CoM_TAU (ndarray): Center of Mass of Pathology Regions (matched for regions we excluded based on observation threshold for TAU)
        pathNames_TAU (str list): list of pathology regions we map to 3D Atlas (matched for regions we excluded based on observation threshold for TAU)
        MarkerVecTAU (ndarray): vector denoting the node size of our 3D Mapped Network (for TAU)
        colorVecTAU (ndarray): vector denoting the node color of our 3D Mapped Network (for TAU)
        CoM_TDP (ndarray): Center of Mass of Pathology Regions (matched for regions we excluded based on observation threshold for TDP)
        pathNames_TDP (str list): list of pathology regions we map to 3D Atlas (matched for regions we excluded based on observation threshold for TDP)
        MarkerVecTDP (ndarray): vector denoting the node size of our 3D Mapped Network (for TDP)
        colorVecTDP (ndarray): vector denoting the node color of our 3D Mapped Network (for TDP)
        cRange (int list): list containing the max min value of Network edge
        outputDir (str): Path location to save the analysis results
        suffix_M (str): Suffix denoting if the analysis is on GM or WM
    """
    # Images to Crop
    TAU_img = None
    TDP_img = None

    TAU_GT_TDP_img = None
    TDP_GT_TAU_img = None

    # Fixed Density Part
    TAU_GT_TDP_FD_img = None
    TDP_GT_TAU_FD_img = None

    # MODIFY the last 3 covMat in covMatlist for FixedDensity
    for i in range(4, 6):
        covMatlist[i] = fixedDensity(covMatlist[i], fd_val)

    covMatNamelist = ['CovMat_FTD_TAU', 'CovMat_FTD_TDP', 'CovMat_FTD_TAU_GT_TDP', 'CovMat_FTD_TDP_GT_TAU', f'CovMat_TAU_gt_TDP_FD_{fd_val}', f'CovMat_TDP_gt_TAU_FD_{fd_val}']
    imglist = [TAU_img, TDP_img, TAU_GT_TDP_img, TDP_GT_TAU_img, TAU_GT_TDP_FD_img, TDP_GT_TAU_FD_img]
    for j in range(len(covMatlist)):
        # Define figure
        fig_atlas = plt.figure()

        if j == 0 or j == 2 or j == 4: # TAU part
            currPathCoM = CoM_TAU
            currLabelNames = pathNames_TAU
            currMarkerVec = MarkerVecTAU
            currColorVec = colorVecTAU
        else: # TDP part
            currPathCoM = CoM_TDP
            currLabelNames = pathNames_TDP
            currMarkerVec = MarkerVecTDP
            currColorVec = colorVecTDP
        
        # Set edge color type (Original - jet / sig - all red)
        if j == 0 or j == 1:
            covType = 'original'
        else:
            covType = 'sig'

        # covMat_TAU, covMat_TDP, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU
        plotNetwork3Individual(fig_atlas, covMatlist[j], NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currPathCoM, 
                               currLabelNames, cRange, currMarkerVec, currColorVec, 3, showLabels = 0, covType = covType)

        fig_atlas.tight_layout()

        # Save Graph Object
        save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + suffix_M)

        # Save the figure
        plt.savefig(outputDir + '/3D_Atlas_' + covMatNamelist[j] + suffix_M, dpi=1000, bbox_inches='tight')

        # Read figure to crop
        imglist[j] = Image.open(outputDir + '/3D_Atlas_' + covMatNamelist[j] + suffix_M + '.png')
        
        
    comb_img = cropImg6(imglist)
    comb_img.save(outputDir + f'/FTD{suffix_M}.png')
####################################################################################################################################################



################################################################## THICKNESS PART ##################################################################
def thickSelectNetwork(networkNames, ii, pathLUT):
    """Select Specific Thickness Network we want to analyze.
    From networkNames=['DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All']

    Args:
        networkNames (str list): list of thickness Network we want to analyze
        ii (int): index to choose from networkNames
        pathLUT (DataFrame): Look up table matching Atlas region names to Atlas Labels (Index) in NetworkDataGeneral Object

    Returns:
        currNetwork (bool list): List of Boolean denoting if Label_Name2_list contains specified networkNames or not / Boolean list of length 400 / in order of Label_Name2_List
        saveName (str): string name of the thickness network we choose
        regNames (str list): Thickness Region Names that are in the specific network we choose (in order of the 400 regions)
        pathNames (str list): Pathology Region Names that matches to regNames (in order of the 400 regions)
        N (int):  Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
    """

    ########################################################### Variable Summary ######################################################################
    # Label_Name2_list: Names of 400 Regions in the Schaffer Atlas (in order with numlab = 400)
    # currNetwork: List of Boolean denoting if Label_Name2_list contains specified networkNames or not / Boolean list of length 400 / in order of Label_Name2_List
    # pathLUTSorted: Sort the pathology LUT by Label_ID400x7 1 to 400 / LUT matching 400 thickness regions to PathName
    # regNames: Thickness Region Names that are in the specific network we choose (in order of the 400 regions)
    # pathNames: Pathology Region Names that matches to regNames (in order of the 400 regions)
    # N: Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
    ###################################################################################################################################################

    # Due to mat file issues, loading Label_Name2 from txt file / list of length 400 --> Names of 400 Regions in the Schaffer Atlas (in order)
    with open('/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/NetworkDataGeneral_Schaefer400x7_LUT_Label_Name2.txt') as file:
        Label_Name2_list = [line.rstrip() for line in file]
        file.close()
    
    if(ii == 8): # Since it's All we set all values to 1 (True) in currNetwork
        currNetwork = np.ones(len(Label_Name2_list), dtype=bool)
    else: # NetworkName --> 'Default' (ii = 5)
        # currNetwork: List of Boolean denoting if Label_Name2_list contains networkNames('Default') or not / Boolean list of length 400
        currNetwork = np.array([networkNames[ii - 1] in label_name for label_name in Label_Name2_list], dtype=bool) # -1 added because of index difference between python and matlab

    # saveName --> this specific case: 'Default' / 'All'
    saveName = networkNames[ii - 1] # -1 added because of index difference between python and matlab

    # Sort the pathology LUT by Label_ID400x7
    pathLUTSorted = pathLUT.sort_values(by=['Label_ID400x7']) # Sort by label 1 to 400
    
    # Get the PathNames of the current Network inside the pathLUTSorted dataset as list / Out of 400 Regions there are 91 regions that contains 'Default'
    pathNames = pathLUTSorted[currNetwork]['PathName'].tolist()
    
    # Get the Label_Name2 that contains current Network ('Default') as list
    regNamesRaw = list(compress(Label_Name2_list, currNetwork))
    
    # N number of Thickness Regions we are analyzing for this specific network
    N = len(regNamesRaw)

    # Change the regNamesRaw (_ to -)
    regNames = []

    for j in range(N):
        nameSplit = regNamesRaw[j].replace('_', '-')
        regNames.append(nameSplit)

    # Sanity Check
    if ii == 5:
        assert(N == 91)
    else: # ii == 8
        assert(N == 400)

    return currNetwork, saveName, regNames, pathNames, N

def generateZScore(HCData, TAUData, TDPData):
    """Generate Z Scores for HC/TAU/TDP Data

    Args:
        HCData (ndarray): HC data
        TAUData (ndarray): TAU data
        TDPData (ndarray): TDP data

    Returns:
        HC_z (float): HC Z Score
        TAU_z (float): TAU Z Score
        TDP_z (float): TDP Z Score
    """
    # Get Mean/STD of the thickness values for HC / Mean, STD of HC for each 400 regions / 1 x 400
    hc_Mean = np.nanmean(HCData, axis=0)
    hc_SD = np.nanstd(HCData, axis=0, ddof=0) # ddof parameter is set to 0, which specifies that the divisor should be N, where N is the number of non-NaN elements in the array

    # HC_z: Z score of Thickness values for Control Dataset
    HC_z = np.empty(HCData.shape)
    for i in range(HCData.shape[0]):
        HC_z[i, :] = (HCData[i, :] - hc_Mean) / hc_SD

    # TAU_z: Z score of Thickness values for Patient Dataset with TAU 
    TAU_z = np.empty(TAUData.shape)
    for i in range(TAUData.shape[0]):
        TAU_z[i, :] = (TAUData[i, :] - hc_Mean) / hc_SD

    # TDP_z: Z score of Thickness values for Patient Dataset with TDP
    TDP_z = np.empty(TDPData.shape)
    for i in range(TDPData.shape[0]):
        TDP_z[i, :] = (TDPData[i, :] - hc_Mean) / hc_SD

    return HC_z, TAU_z, TDP_z

def generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCData, TAUData, TDPData):
    """Generate W Scores for HC/TAU/TDP Data

    Args:
        AgeSexHC (ndarray): Age and Sex data of HC 
        AgeSexTAU (ndarray): Age and Sex data of TAU 
        AgeSexTDP (ndarray): Age and Sex data of TDP 
        N (int):  Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
        HCData (ndarray): HC data
        TAUData (ndarray): TAU data
        TDPData (ndarray): TDP data

    Returns:
        HC_w (float): HC W Score
        TAU_w (float): TAU W Score
        TDP_w (float): TDP W Score
    """

    HC_model_list = [] # (400,) 
    HC_residuals_std_list = [] # (400,)  

    # Linear regression on HC for ALL 400 Regions
    for k in range(N): 
        yHC = HCData[:,k] # Thickness values of specified region for HC (54,)
        assert(yHC.shape == (HCData.shape[0],))
    
        # Linear Regression
        regHC = LinearRegression().fit(AgeSexHC, yHC)

        # Predict 
        HC_Predict = regHC.predict(AgeSexHC) # Shape (54, ) --> Thickness values for 54 subject in specified Region
        assert(HC_Predict.shape == (HCData.shape[0],))

        # Residual
        HC_Residual = yHC - HC_Predict # Shape (54,)
        assert(HC_Residual.shape == (HCData.shape[0],))

        HC_std = np.nanstd(HC_Residual, axis=0, ddof=0)

        # Save the Linear Regression Model
        HC_model_list.append(regHC)

        # Save residual std
        HC_residuals_std_list.append(HC_std)

    # Sanity Check
    assert(len(HC_model_list) == N) # list length of 400
    assert(len(HC_residuals_std_list) == N) # list length of 400

    # Predict values of HC, TAU and TDP using these coef / (54, 400), (26, 400), (30, 400)
    HC_Predicted = np.empty(HCData.shape) # (54, 400) 
    for h in range(HCData.shape[0]): # 54
        for g in range (HCData.shape[1]): # 400 or 40
            HC_Predicted[h, g] = HC_model_list[g].predict([AgeSexHC[h]])

    TAU_Predicted = np.empty(TAUData.shape) # (26, 400)
    for h in range(TAUData.shape[0]): # 26
        for g in range (TAUData.shape[1]): # 400 or 35
            TAU_Predicted[h, g] = HC_model_list[g].predict([AgeSexTAU[h]])

    TDP_Predicted = np.empty(TDPData.shape) # (30, 400) 
    for h in range(TDPData.shape[0]): # 30
        for g in range (TDPData.shape[1]): # 400 or 35
            TDP_Predicted[h, g] = HC_model_list[g].predict([AgeSexTDP[h]])

    # Compute W Score
    HC_w = np.empty(HCData.shape) # 54 x 400
    TAU_w = np.empty(TAUData.shape) # 26 x 400
    TDP_w = np.empty(TDPData.shape) # 30 x 400

    for i in range(HCData.shape[0]): # 54
        HC_w[i, :] = (HCData[i, :] - HC_Predicted[i, :]) / HC_residuals_std_list

    for i in range(TAUData.shape[0]): # 26
        TAU_w[i, :] = (TAUData[i, :] - TAU_Predicted[i, :]) / HC_residuals_std_list

    for i in range(TDPData.shape[0]): # 30
        TDP_w[i, :] = (TDPData[i, :] - TDP_Predicted[i, :]) / HC_residuals_std_list

    return HC_w, TAU_w, TDP_w

def sexLUT(sexstr):
    """convert string (gender type) to int

    Args:
        sexstr (str): Male or Female

    Raises:
        Exception: Other than Male or Female

    Returns:
        0 or 1 (int): 0 - Male / 1 - Female
    """
    if sexstr == 'Male':
        return 0
    elif sexstr == 'Female':
        return 1
    else:
        raise Exception("Problem with converting sex to integer")

def sexToBinary(listOfValues):
    return [sexLUT(i) for i in listOfValues]

def thicknessCovMatGenerate_SigEdgesTSave(N, HCData, TAUData, TDPData, pthresh, cov_thresh, outputDir, networkNames, regNames, pathNames, ii):
    """Calculate cov Mat for thickness data analysis (total of 15 cov Mat) - with additional saving SigEdgesT.csv

    Args:
        N (int):  Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
        HCData (ndarray): HC data
        TAUData (ndarray): TAU data
        TDPData (ndarray): TDP data
        pthresh (float): p-value threshold we are currently using
        cov_thresh (float): covariance threshold we are currently using (to remove potential noise)
        outputDir (str): Path location to save the analysis results
        networkNames (str list): list of thickness Network we want to analyze
        regNames (str list): Thickness Region Names that are in the specific network we choose (in order of the 400 regions)
        pathNames (str list): Pathology Region Names that matches to regNames (in order of the 400 regions)
        ii (int): index to choose from networkNames
    Returns:
        (list of ndarray): return list of 15 cov Mat calculated
    """
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

    # N --> Number of Regions to generate cov Mat
    covMatHC = np.full((N,N), np.nan)
    covMatTAU = np.full((N,N), np.nan)
    covMatTDP = np.full((N,N), np.nan)

    cmpCovTAU_gt_HC = np.full((N,N), np.nan)
    cmpCovTDP_gt_HC = np.full((N,N), np.nan)

    cmpCovTAU_lt_HC = np.full((N,N), np.nan)
    cmpCovTDP_lt_HC = np.full((N,N), np.nan)

    cmpCovTAU_gt_TDP = np.full((N,N), np.nan)
    cmpCovTDP_gt_TAU = np.full((N,N), np.nan)

    cmpCovTAU_gt_HC_raw = np.full((N,N), np.nan)
    cmpCovTDP_gt_HC_raw = np.full((N,N), np.nan)

    cmpCovTAU_lt_HC_raw = np.full((N,N), np.nan)
    cmpCovTDP_lt_HC_raw = np.full((N,N), np.nan)

    cmpCovTAU_gt_TDP_raw = np.full((N,N), np.nan)
    cmpCovTDP_gt_TAU_raw = np.full((N,N), np.nan)

    for i in range(N):
        for j in range(N):

            if i != j: # When it's not the same region in the current Network
                nHC = np.sum(~np.isnan(HCData[:, i]) & ~np.isnan(HCData[:, j])) # Get the sum of regions that are not NaN
                nTAU = np.sum(~np.isnan(TAUData[:, i]) & ~np.isnan(TAUData[:, j])) # Get the sum of regions that are not NaN
                nTDP = np.sum(~np.isnan(TDPData[:, i]) & ~np.isnan(TDPData[:, j])) # Get the sum of regions that are not NaN

                if nHC > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatHC[i, j] = ma.corrcoef(ma.masked_invalid(HCData[:,i]), ma.masked_invalid(HCData[:,j]))[0, 1]
                    if covMatHC[i, j] < cov_thresh:
                        covMatHC[i, j] = np.nan
                        nHC = 0 # reset for below case

                if nTAU > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTAU[i, j] = ma.corrcoef(ma.masked_invalid(TAUData[:,i]), ma.masked_invalid(TAUData[:,j]))[0, 1]
                    if covMatTAU[i, j] < cov_thresh:
                        covMatTAU[i, j] = np.nan
                        nTAU = 0

                if nTDP > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covMatTDP[i, j] < cov_thresh:
                        covMatTDP[i, j] = np.nan
                        nTDP = 0

                if nHC > 3 and nTAU > 3: 
                    cmpCovTAU_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTAU[i, j], nHC, nTAU)[0] < pthresh # get only probs
                    cmpCovTAU_lt_HC[i, j] = test_corr_sig(covMatTAU[i, j], covMatHC[i, j], nTAU, nHC)[0] < pthresh # get only probs

                    cmpCovTAU_gt_HC_raw[i, j] = test_corr_sig(covMatHC[i, j], covMatTAU[i, j], nHC, nTAU)[0] # get raw sig value
                    cmpCovTAU_lt_HC_raw[i, j] = test_corr_sig(covMatTAU[i, j], covMatHC[i, j], nTAU, nHC)[0] # get raw sig value

                if nHC > 3 and nTDP > 3: 
                    cmpCovTDP_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTDP[i, j], nHC, nTDP)[0] < pthresh # get only probs
                    cmpCovTDP_lt_HC[i, j] = test_corr_sig(covMatTDP[i, j], covMatHC[i, j], nTDP, nHC)[0] < pthresh # get only probs

                    cmpCovTDP_gt_HC_raw[i, j] = test_corr_sig(covMatHC[i, j], covMatTDP[i, j], nHC, nTDP)[0] 
                    cmpCovTDP_lt_HC_raw[i, j] = test_corr_sig(covMatTDP[i, j], covMatHC[i, j], nTDP, nHC)[0] 

                if nTAU > 3 and nTDP > 3:
                    cmpCovTAU_gt_TDP[i, j] = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)[0] < pthresh # get only probs
                    cmpCovTDP_gt_TAU[i, j] = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)[0] < pthresh # get only probs

                    cmpCovTAU_gt_TDP_raw[i, j] = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)[0] 
                    cmpCovTDP_gt_TAU_raw[i, j] = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)[0] 
                

                # Fill in the SigEdgesT pandas dataframe we created earlier
                # When there is a significance diff between TAU AND TDP
                if (cmpCovTAU_gt_TDP[i, j] == 1 or cmpCovTDP_gt_TAU[i, j] == 1) and i < j: # record significant edges
                    SigEdgesCounter += 1 # increase SigEdgesCounter
                    ss = SigEdgesCounter - 1 # index difference between matlab and python

                    SigEdgesT.at[ss, 'NetworkNames'] = networkNames[ii - 1] # index difference between matlab and python
                    SigEdgesT.at[ss, 'ROIName1'] = regNames[i]
                    SigEdgesT.at[ss, 'ROIName1_Path'] = pathNames[i] # NaN is obmitted         
                    SigEdgesT.at[ss, 'ROIName2'] = regNames[j]
                    SigEdgesT.at[ss, 'ROIName2_Path'] = pathNames[j] # NaN is obmitted 
                    SigEdgesT.at[ss, 'r_YesAD'] = ma.corrcoef(ma.masked_invalid(TAUData[:,i]), ma.masked_invalid(TAUData[:,j]))[0, 1]
                    SigEdgesT.at[ss, 'r_NoAD'] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    SigEdgesT.at[ss, 'N_YesAD'] = nTAU
                    SigEdgesT.at[ss, 'N_NoAD'] = nTDP

                    p, z = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)
                    SigEdgesT.at[ss, 'p_YesAD_gt_NoAD'] = p
                    SigEdgesT.at[ss, 'z_YesAD_gt_NoAD'] = z

                    p, z = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)
                    SigEdgesT.at[ss, 'p_NoAD_gt_YesAD'] = p
                    SigEdgesT.at[ss, 'z_NoAD_gt_YesAD'] = z
    

    # Save SigEdgesT as csv
    SigEdgesT.to_csv(outputDir + "/SigEdges_FTD.csv")  

    return [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU, cmpCovTAU_gt_HC_raw, cmpCovTDP_gt_HC_raw, cmpCovTAU_lt_HC_raw, cmpCovTDP_lt_HC_raw, cmpCovTAU_gt_TDP_raw, cmpCovTDP_gt_TAU_raw]

def thicknessCovMatGenerate(N, HCData, TAUData, TDPData, pthresh, cov_thresh):
    """Calculate cov Mat for thickness data analysis (total of 15 cov Mat)

    Args:
        N (int):  Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
        HCData (ndarray): HC data
        TAUData (ndarray): TAU data
        TDPData (ndarray): TDP data
        pthresh (float): p-value threshold we are currently using
        cov_thresh (float): covariance threshold we are currently using (to remove potential noise)

    Returns:
        (list of ndarray): return list of 15 cov Mat calculated
    """
    # N --> Number of Regions to generate cov Mat
    covMatHC = np.full((N,N), np.nan)
    covMatTAU = np.full((N,N), np.nan)
    covMatTDP = np.full((N,N), np.nan)

    cmpCovTAU_gt_HC = np.full((N,N), np.nan)
    cmpCovTDP_gt_HC = np.full((N,N), np.nan)

    cmpCovTAU_lt_HC = np.full((N,N), np.nan)
    cmpCovTDP_lt_HC = np.full((N,N), np.nan)

    cmpCovTAU_gt_TDP = np.full((N,N), np.nan)
    cmpCovTDP_gt_TAU = np.full((N,N), np.nan)

    cmpCovTAU_gt_HC_raw = np.full((N,N), np.nan)
    cmpCovTDP_gt_HC_raw = np.full((N,N), np.nan)

    cmpCovTAU_lt_HC_raw = np.full((N,N), np.nan)
    cmpCovTDP_lt_HC_raw = np.full((N,N), np.nan)

    cmpCovTAU_gt_TDP_raw = np.full((N,N), np.nan)
    cmpCovTDP_gt_TAU_raw = np.full((N,N), np.nan)

    for i in range(N):
        for j in range(N):

            if i != j: # When it's not the same region in the current Network
                nHC = np.sum(~np.isnan(HCData[:, i]) & ~np.isnan(HCData[:, j])) # Get the sum of regions that are not NaN
                nTAU = np.sum(~np.isnan(TAUData[:, i]) & ~np.isnan(TAUData[:, j])) # Get the sum of regions that are not NaN
                nTDP = np.sum(~np.isnan(TDPData[:, i]) & ~np.isnan(TDPData[:, j])) # Get the sum of regions that are not NaN

                if nHC > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatHC[i, j] = ma.corrcoef(ma.masked_invalid(HCData[:,i]), ma.masked_invalid(HCData[:,j]))[0, 1]
                    if covMatHC[i, j] < cov_thresh:
                        covMatHC[i, j] = np.nan
                        nHC = 0 # reset for below case

                if nTAU > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTAU[i, j] = ma.corrcoef(ma.masked_invalid(TAUData[:,i]), ma.masked_invalid(TAUData[:,j]))[0, 1]
                    if covMatTAU[i, j] < cov_thresh:
                        covMatTAU[i, j] = np.nan
                        nTAU = 0 

                if nTDP > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covMatTDP[i, j] < cov_thresh:
                        covMatTDP[i, j] = np.nan
                        nTDP = 0

                if nHC > 3 and nTAU > 3: 
                    cmpCovTAU_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTAU[i, j], nHC, nTAU)[0] < pthresh # get only probs
                    cmpCovTAU_lt_HC[i, j] = test_corr_sig(covMatTAU[i, j], covMatHC[i, j], nTAU, nHC)[0] < pthresh # get only probs

                    cmpCovTAU_gt_HC_raw[i, j] = test_corr_sig(covMatHC[i, j], covMatTAU[i, j], nHC, nTAU)[0] # get raw sig value
                    cmpCovTAU_lt_HC_raw[i, j] = test_corr_sig(covMatTAU[i, j], covMatHC[i, j], nTAU, nHC)[0] # get raw sig value

                if nHC > 3 and nTDP > 3: 
                    cmpCovTDP_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTDP[i, j], nHC, nTDP)[0] < pthresh # get only probs
                    cmpCovTDP_lt_HC[i, j] = test_corr_sig(covMatTDP[i, j], covMatHC[i, j], nTDP, nHC)[0] < pthresh # get only probs

                    cmpCovTDP_gt_HC_raw[i, j] = test_corr_sig(covMatHC[i, j], covMatTDP[i, j], nHC, nTDP)[0] 
                    cmpCovTDP_lt_HC_raw[i, j] = test_corr_sig(covMatTDP[i, j], covMatHC[i, j], nTDP, nHC)[0] 

                if nTAU > 3 and nTDP > 3:
                    cmpCovTAU_gt_TDP[i, j] = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)[0] < pthresh # get only probs
                    cmpCovTDP_gt_TAU[i, j] = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)[0] < pthresh # get only probs

                    cmpCovTAU_gt_TDP_raw[i, j] = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)[0] 
                    cmpCovTDP_gt_TAU_raw[i, j] = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)[0] 

    return [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU, cmpCovTAU_gt_HC_raw, cmpCovTDP_gt_HC_raw, cmpCovTAU_lt_HC_raw, cmpCovTDP_lt_HC_raw, cmpCovTAU_gt_TDP_raw, cmpCovTDP_gt_TAU_raw]

def fixedDensity(covMat, N):

    # Find the smallest Nth value. 
    flatArray = covMat.flatten()

    # Sort the array
    flatArray_Sorted = np.sort(flatArray)

    # Get the Nth smallest value (P-value)
    N = 2 * N # Because it is adjacency matrix
    threshVal = flatArray_Sorted[N-1]

    # Return a CovMat, that only keeps the smallest N values. Everything else --> NaN
    fixedDcovMat = np.where(covMat <= threshVal, covMat, np.nan) # True / False

    return fixedDcovMat

def thickness3DMapping(covMatlist, NetworkDataGeneral, currCoM, LabelNames, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVec, outputDir, pthresh, fd_val, FD = True, W_Score = False):
    """3D mapping of the Pathology Data Analysis

    Args:
        covMatlist (list of ndarray): list of 15 cov mat we calculated for thickness dataset
        NetworkDataGeneral (obj): Preconstructed atlas data
        currCoM ndarray): Center of Mass of Thickness Regions
        LabelNames (str list): list of thickness regions we map to 3D Atlas
        cRange (int list): list containing the max min value of Network edge
        MarkerVecHC (ndarray): vector denoting the node size of our 3D Mapped Network (for HC)
        MarkerVecTAU (ndarray): vector denoting the node size of our 3D Mapped Network (for TAU)
        MarkerVecTDP (ndarray): vector denoting the node size of our 3D Mapped Network (for TDP)
        colorVec (ndarray): vector denoting the node color of our 3D Mapped Network (for HC, TAU, TDP)
        outputDir (str): Path location to save the analysis results
        pthresh (float): p-value threshold we are currently using
        fd_val (int): Fixed Density value. 100 --> we only look at top 100 most significant edges (lowest p-values)
        FD (bool, optional): boolean denoting if we want to do FD analysis. Defaults to True.
        W_Score (bool, optional): boolean denoting if this is W-Score 3D Mapping or not. Defaults to False.
    """
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
    # Fixed Density Part
    TAU_gt_HC_FD_img = None
    TDP_gt_HC_FD_img = None
    TAU_lt_HC_FD_img = None
    TDP_lt_HC_FD_img = None
    TAU_gt_TDP_FD_img = None
    TDP_gt_TAU_FD_img = None
    
    imglist = [HC_img, TAU_img, TDP_img, TAU_gt_HC_img, TDP_gt_HC_img, TAU_lt_HC_img, TDP_lt_HC_img, TAU_gt_TDP_img, TDP_gt_TAU_img, TAU_gt_HC_FD_img, TDP_gt_HC_FD_img, TAU_lt_HC_FD_img, TDP_lt_HC_FD_img, TAU_gt_TDP_FD_img, TDP_gt_TAU_FD_img]

    # MODIFY the last 6 covMat in covMatlist for FixedDensity
    for i in range(9, 15):
        covMatlist[i] = fixedDensity(covMatlist[i], fd_val)
    
    if W_Score: # W_Score
        covMatNamelist = ['covMatHC_W', 'covMatTAU_W', 'covMatTDP_W', 'cmpCovTAU_gt_HC_W', 'cmpCovTDP_gt_HC_W', 'cmpCovTAU_lt_HC_W', 'cmpCovTDP_lt_HC_W', 'cmpCovTAU_gt_TDP_W', 'cmpCovTDP_gt_TAU_W', f'cmpCovTAU_gt_HC_FD_W_{fd_val}', f'cmpCovTDP_gt_HC_FD_W_{fd_val}', f'cmpCovTAU_lt_HC_FD_W_{fd_val}', f'cmpCovTDP_lt_HC_FD_W_{fd_val}', f'cmpCovTAU_gt_TDP_FD_W_{fd_val}', f'cmpCovTDP_gt_TAU_FD_W_{fd_val}']
    else : # Original
        covMatNamelist = ['covMatHC', 'covMatTAU', 'covMatTDP', 'cmpCovTAU_gt_HC', 'cmpCovTDP_gt_HC', 'cmpCovTAU_lt_HC', 'cmpCovTDP_lt_HC', 'cmpCovTAU_gt_TDP', 'cmpCovTDP_gt_TAU', f'cmpCovTAU_gt_HC_FD_{fd_val}', f'cmpCovTDP_gt_HC_FD_{fd_val}', f'cmpCovTAU_lt_HC_FD_{fd_val}', f'cmpCovTDP_lt_HC_FD_{fd_val}', f'cmpCovTAU_gt_TDP_FD_{fd_val}', f'cmpCovTDP_gt_TAU_FD_{fd_val}']
    
        
    MarkerVecList = [MarkerVecHC, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP, MarkerVecTAU, MarkerVecTDP]
    
    for j in range(len(covMatlist)):
        # Define figure
        fig_atlas = plt.figure()

        # Define MarkerVec to Use / Same np.ones(np.sum(currNetwork))
        MarkerVec = MarkerVecList[j]

        # Edge color colormap
        if j == 0 or j == 1 or j == 2:
            covType = 'original'
        else:
            covType = 'sig'

        # 3D Mapping
        plotNetwork3Individual(fig_atlas, covMatlist[j], NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MarkerVec, colorVec, 3, showLabels = 0, covType = covType)

        fig_atlas.tight_layout()

        if W_Score: # W_score
            # Save Graph Object
            save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_WScore')
            # Save the figure
            plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}', dpi=1000, bbox_inches='tight')
            # Read figure to crop
            imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}' + '.png')
        else: # Original
            # Save Graph Object
            save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}')
            # Save the figure
            plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}', dpi=1000, bbox_inches='tight')
            # Read figure to crop
            imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}' + '.png')


    thickness_comb_img = cropImg6(imglist[3:9])
    thickness_comb_FD_img = cropImg6(imglist[9:])

    if W_Score: # W_Score
        thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(WSCORE_Original).png')
        thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(WSCORE_FixedDensity).png')
    else: # Original
        thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(Original).png')
        thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}(FixedDensity).png')
####################################################################################################################################################

############################################################### THICKNESS AT PATH PART ##############################################################
def thicknessAtPathFewObs(covMatlist, TAU_missing_index, TDP_missing_index):
    # Modify the covMatlist data so exclude these regions
    
    # Rows
    # covMatlist[0] = np.delete(covMatlist[0], TAU_missing_index, axis = 0) # covMatHC
    covMatlist[1] = np.delete(covMatlist[1], TAU_missing_index, axis = 0) # covMatTAU
    covMatlist[2] = np.delete(covMatlist[2], TDP_missing_index, axis = 0) # covMatTDP
    covMatlist[3] = np.delete(covMatlist[3], TAU_missing_index, axis = 0) # cmpCovTAU_gt_HC
    covMatlist[4] = np.delete(covMatlist[4], TDP_missing_index, axis = 0) # cmpCovTDP_gt_HC
    covMatlist[5] = np.delete(covMatlist[5], TAU_missing_index, axis = 0) # cmpCovTAU_lt_HC
    covMatlist[6] = np.delete(covMatlist[6], TDP_missing_index, axis = 0) # cmpCovTDP_lt_HC
    covMatlist[7] = np.delete(covMatlist[7], TAU_missing_index, axis = 0) # cmpCovTAU_gt_TDP
    covMatlist[8] = np.delete(covMatlist[8], TDP_missing_index, axis = 0) # cmpCovTDP_gt_TAU
    covMatlist[9] = np.delete(covMatlist[9], TAU_missing_index, axis = 0) # cmpCovTAU_gt_HC_raw
    covMatlist[10] = np.delete(covMatlist[10], TDP_missing_index, axis = 0) # cmpCovTDP_gt_HC_raw
    covMatlist[11] = np.delete(covMatlist[11], TAU_missing_index, axis = 0) # cmpCovTAU_lt_HC_raw
    covMatlist[12] = np.delete(covMatlist[12], TDP_missing_index, axis = 0) # cmpCovTDP_lt_HC_raw
    covMatlist[13] = np.delete(covMatlist[13], TAU_missing_index, axis = 0) # cmpCovTAU_gt_TDP_raw
    covMatlist[14] = np.delete(covMatlist[14], TDP_missing_index, axis = 0) # cmpCovTDP_gt_TAU_raw

    # Columns
    # covMatlist[0] = np.delete(covMatlist[0], TAU_missing_index, axis = 0) # covMatHC
    covMatlist[1] = np.delete(covMatlist[1], TAU_missing_index, axis = 1) # covMatTAU
    covMatlist[2] = np.delete(covMatlist[2], TDP_missing_index, axis = 1) # covMatTDP
    covMatlist[3] = np.delete(covMatlist[3], TAU_missing_index, axis = 1) # cmpCovTAU_gt_HC
    covMatlist[4] = np.delete(covMatlist[4], TDP_missing_index, axis = 1) # cmpCovTDP_gt_HC
    covMatlist[5] = np.delete(covMatlist[5], TAU_missing_index, axis = 1) # cmpCovTAU_lt_HC
    covMatlist[6] = np.delete(covMatlist[6], TDP_missing_index, axis = 1) # cmpCovTDP_lt_HC
    covMatlist[7] = np.delete(covMatlist[7], TAU_missing_index, axis = 1) # cmpCovTAU_gt_TDP
    covMatlist[8] = np.delete(covMatlist[8], TDP_missing_index, axis = 1) # cmpCovTDP_gt_TAU
    covMatlist[9] = np.delete(covMatlist[9], TAU_missing_index, axis = 1) # cmpCovTAU_gt_HC_raw
    covMatlist[10] = np.delete(covMatlist[10], TDP_missing_index, axis = 1) # cmpCovTDP_gt_HC_raw
    covMatlist[11] = np.delete(covMatlist[11], TAU_missing_index, axis = 1) # cmpCovTAU_lt_HC_raw
    covMatlist[12] = np.delete(covMatlist[12], TDP_missing_index, axis = 1) # cmpCovTDP_lt_HC_raw
    covMatlist[13] = np.delete(covMatlist[13], TAU_missing_index, axis = 1) # cmpCovTAU_gt_TDP_raw
    covMatlist[14] = np.delete(covMatlist[14], TDP_missing_index, axis = 1) # cmpCovTDP_gt_TAU_raw

    return covMatlist

def thicknessAtPath3DMapping(covMatlist, NetworkDataGeneral, CoM_HC, CoM_TAU, CoM_TDP, LabelNames_Path_HC, LabelNames_Path_TAU, LabelNames_Path_TDP, cRange, MarkerVecHC, MarkerVecTAU, MarkerVecTDP, colorVecHC, colorVecTAU, colorVecTDP, outputDir, pthresh, fd_val, FD = True, W_Score = False):
    
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

    # MODIFY the last 6 covMat in covMatlist for FixedDensity
    for i in range(9, 15):
        covMatlist[i] = fixedDensity(covMatlist[i], fd_val)

    if W_Score: # W_Score
        covMatNamelist = ['covMatHC_Path_W', 'covMatTAU_Path_W', 'covMatTDP_Path_W', 'cmpCovTAU_gt_HC_Path_W', 'cmpCovTDP_gt_HC_Path_W', 'cmpCovTAU_lt_HC_Path_W', 'cmpCovTDP_lt_HC_Path_W', 'cmpCovTAU_gt_TDP_Path_W', 'cmpCovTDP_gt_TAU_Path_W', f'cmpCovTAU_gt_HC_FD_Path_W_{fd_val}', f'cmpCovTDP_gt_HC_FD_Path_W_{fd_val}', f'cmpCovTAU_lt_HC_FD_Path_W_{fd_val}', f'cmpCovTDP_lt_HC_FD_Path_W_{fd_val}', f'cmpCovTAU_gt_TDP_FD_Path_W_{fd_val}', f'cmpCovTDP_gt_TAU_FD_Path_W_{fd_val}']
    else: # Original
        covMatNamelist = ['covMatHC_Path', 'covMatTAU_Path', 'covMatTDP_Path', 'cmpCovTAU_gt_HC_Path', 'cmpCovTDP_gt_HC_Path', 'cmpCovTAU_lt_HC_Path', 'cmpCovTDP_lt_HC_Path', 'cmpCovTAU_gt_TDP_Path', 'cmpCovTDP_gt_TAU_Path', f'cmpCovTAU_gt_HC_FD_{fd_val}_Path', f'cmpCovTDP_gt_HC_FD_{fd_val}_Path', f'cmpCovTAU_lt_HC_FD_{fd_val}_Path', f'cmpCovTDP_lt_HC_FD_{fd_val}_Path', f'cmpCovTAU_gt_TDP_FD_{fd_val}_Path', f'cmpCovTDP_gt_TAU_FD_{fd_val}_Path']

    for j in range(len(covMatlist)):
        # Define figure
        fig_atlas = plt.figure()

        # Edge color colormap
        if j == 0 or j == 1 or j == 2:
            covType = 'original'
        else:
            covType = 'sig'

        if j == 0: # HC
            currPathCoM = CoM_HC
            currLabelNames_Path = LabelNames_Path_HC
            currMarkerVec = MarkerVecHC
            currColorVec = colorVecHC
        elif j == 1 or j == 3 or j == 5 or j == 7: # TAU part
            currPathCoM = CoM_TAU
            currLabelNames_Path = LabelNames_Path_TAU
            currMarkerVec = MarkerVecTAU
            currColorVec = colorVecTAU
        else: # TDP part
            currPathCoM = CoM_TDP
            currLabelNames_Path = LabelNames_Path_TDP
            currMarkerVec = MarkerVecTDP
            currColorVec = colorVecTDP

        # [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]
        plotNetwork3Individual(fig_atlas, covMatlist[j], NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currPathCoM, currLabelNames_Path, cRange, currMarkerVec, currColorVec, 3, showLabels = 0, covType = covType)

        fig_atlas.tight_layout()

        if W_Score: # W_Score
            # Save Graph Object
            save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_Path_pthresh_{str(pthresh).split(".")[1]}_WScore')
            # Save the figure
            plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)', dpi=1000, bbox_inches='tight')
            # Read figure to crop
            imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + '_WSCORE' + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)' + '.png')

        else: # Original
            # Save Graph Object
            save3DGraph(covMatlist[j], outputDir, covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}')
            # Save the figure
            plt.savefig(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)', dpi=1000, bbox_inches='tight')
            # Read figure to crop
            imglist[j] = Image.open(outputDir + '/Thickness_3D_Atlas_' + covMatNamelist[j] + f'_pthresh_{str(pthresh).split(".")[1]}_(At Path)' + '.png')

    thickness_comb_img = cropImg6(imglist[3:9])
    thickness_comb_FD_img = cropImg6(imglist[9:])

    if W_Score: # W_Score
        thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(WSCORE).png')
        thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(WSCORE_FixedDensity).png')
    else: # Original
        thickness_comb_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(Original).png')
        thickness_comb_FD_img.save(outputDir + f'/Combined_FTD_Thickness_pthresh_{str(pthresh).split(".")[1]}_AT_PATH(FixedDensity).png')

####################################################################################################################################################


def sliceCovMat(covMat, colList):
    return covMat[np.ix_(colList, colList)]

