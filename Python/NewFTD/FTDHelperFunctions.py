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

def sexLUT(sexstr):
    if sexstr == 'Male':
        return 0
    elif sexstr == 'Female':
        return 1
    else:
        raise Exception("Problem with converting sex to integer")

def sexToBinary(listOfValues):
    return [sexLUT(i) for i in listOfValues]

def sliceCovMat(covMat, colList):
    return covMat[np.ix_(colList, colList)]

# Helper Functions
def test_corr_sig(corr1, corr2, ncorr1, ncorr2):
    z_corr1 = np.arctanh(corr1)
    z_corr2 = np.arctanh(corr2)
    z_obs = (z_corr1 - z_corr2) / ((1 / (ncorr1 - 3)) + (1 / (ncorr2 - 3))) ** 0.5
    probs = norm.cdf(z_obs)
    return probs, z_obs

def thicknessCovMatGenerate_SigEdgesTSave(N, HCData, TAUData, TDPData, pthresh, cov_thresh, outputDir, networkNames, regNames, pathNames, ii):
    
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

                if nTDP > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covMatTDP[i, j] < cov_thresh:
                        covMatTDP[i, j] = np.nan

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

                if nTDP > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covMatTDP[i, j] < cov_thresh:
                        covMatTDP[i, j] = np.nan

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

def thicknessCovMatGenerateAtPath(N, HCData, TAUData, TDPData, pthresh, cov_thresh):

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

                if nTDP > 3: # For only where regions that doesn't have invalid values / Nan for where covariance is less than 0.1
                    covMatTDP[i, j] = ma.corrcoef(ma.masked_invalid(TDPData[:,i]), ma.masked_invalid(TDPData[:,j]))[0, 1]
                    if covMatTDP[i, j] < cov_thresh:
                        covMatTDP[i, j] = np.nan

                if nHC > 3 and nTAU > 3: 
                    cmpCovTAU_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTAU[i, j], nHC, nTAU)[0] < pthresh # get only probs
                    cmpCovTAU_lt_HC[i, j] = test_corr_sig(covMatTAU[i, j], covMatHC[i, j], nTAU, nHC)[0] < pthresh # get only probs

                   

                if nHC > 3 and nTDP > 3: 
                    cmpCovTDP_gt_HC[i, j] = test_corr_sig(covMatHC[i, j], covMatTDP[i, j], nHC, nTDP)[0] < pthresh # get only probs
                    cmpCovTDP_lt_HC[i, j] = test_corr_sig(covMatTDP[i, j], covMatHC[i, j], nTDP, nHC)[0] < pthresh # get only probs

           

                if nTAU > 3 and nTDP > 3:
                    cmpCovTAU_gt_TDP[i, j] = test_corr_sig(covMatTDP[i, j], covMatTAU[i, j], nTDP, nTAU)[0] < pthresh # get only probs
                    cmpCovTDP_gt_TAU[i, j] = test_corr_sig(covMatTAU[i, j], covMatTDP[i, j], nTAU, nTDP)[0] < pthresh # get only probs


    return [covMatHC, covMatTAU, covMatTDP, cmpCovTAU_gt_HC, cmpCovTDP_gt_HC, cmpCovTAU_lt_HC, cmpCovTDP_lt_HC, cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU]


def generateZScore(HCData, TAUData, TDPData, hc_Mean, hc_SD):
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

#### NOT USED ####
def generateZScore_AtPath(HCData, TAUData, TDPData, hc_Mean, hc_SD, TAU_missing_index_GM, TDP_missing_index_GM):
    # HC_z: Z score of Thickness values for Control Dataset
    HC_z = np.empty(HCData.shape)
    for i in range(HCData.shape[0]):
        HC_z[i, :] = (HCData[i, :] - hc_Mean) / hc_SD

    # TAU_z: Z score of Thickness values for Patient Dataset with TAU 
    TAU_z = np.empty(TAUData.shape)
    for i in range(TAUData.shape[0]):
        TAU_z[i, :] = (TAUData[i, :] - np.delete(hc_Mean, TAU_missing_index_GM, axis = 0)) / np.delete(hc_SD, TAU_missing_index_GM, axis = 0)

    # TDP_z: Z score of Thickness values for Patient Dataset with TDP
    TDP_z = np.empty(TDPData.shape)
    for i in range(TDPData.shape[0]):
        TDP_z[i, :] = (TDPData[i, :] - np.delete(hc_Mean, TDP_missing_index_GM, axis = 0)) / np.delete(hc_SD, TDP_missing_index_GM, axis = 0)

    return HC_z, TAU_z, TDP_z

def generateWScore_AtPath(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCData, TAUData, TDPData, TAU_missing_index_GM, TDP_missing_index_GM):
    
    HC_model_list = [] # (40,)
    HC_residuals_std_list = [] # (40,)

    # Linear regression on HC for ALL 40 Regions
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

    
    # Convert HC_model_list to numpy array
    HC_model_list = np.array(HC_model_list)
    # Modify HC_model_list to match the Pathology Regions (Removed NaN regions)
    HC_model_list_TAU = np.delete(HC_model_list, TAU_missing_index_GM, axis = 0)
    HC_model_list_TDP = np.delete(HC_model_list, TDP_missing_index_GM, axis = 0)

    # Convert HC_residuals_std_list to numpy array
    HC_residuals_std_list = np.array(HC_residuals_std_list)
    # Modify HC_model_list to match the Pathology Regions (Removed NaN regions)
    HC_residuals_std_list_TAU = np.delete(HC_residuals_std_list, TAU_missing_index_GM, axis = 0)
    HC_residuals_std_list_TDP = np.delete(HC_residuals_std_list, TDP_missing_index_GM, axis = 0)

    # Sanity Check
    assert(HC_model_list.shape == (HCData.shape[1],)) # list length of 40
    assert(HC_model_list_TAU.shape == (TAUData.shape[1],)) # list length of 35
    assert(HC_model_list_TDP.shape == (TDPData.shape[1],)) # list length of 34
    
    assert(HC_residuals_std_list.shape == (HCData.shape[1],)) # list length of 40
    assert(HC_residuals_std_list_TAU.shape == (TAUData.shape[1],)) # list length of 35
    assert(HC_residuals_std_list_TDP.shape == (TDPData.shape[1],)) # list length of 34

    # Predict values of HC, TAU and TDP using these coef / (54, 40), (26, 35), (30, 34)
    HC_Predicted = np.empty(HCData.shape) # (54, 40)
    for h in range(HCData.shape[0]): # 54
        for g in range (HCData.shape[1]): # 40
            HC_Predicted[h, g] = HC_model_list[g].predict([AgeSexHC[h]])

    TAU_Predicted = np.empty(TAUData.shape) # (26, 35)
    for h in range(TAUData.shape[0]): # 26
        for g in range (TAUData.shape[1]): # 35
            TAU_Predicted[h, g] = HC_model_list_TAU[g].predict([AgeSexTAU[h]])

    TDP_Predicted = np.empty(TDPData.shape) # (30, 400) or (30, 34)
    for h in range(TDPData.shape[0]): # 30
        for g in range (TDPData.shape[1]): # 34
            TDP_Predicted[h, g] = HC_model_list_TDP[g].predict([AgeSexTDP[h]])

    # Compute W Score
    HC_w = np.empty(HCData.shape) # 54 x 40
    TAU_w = np.empty(TAUData.shape) # 26 x 35
    TDP_w = np.empty(TDPData.shape) # 30 x 34

    for i in range(HCData.shape[0]): # 54
        HC_w[i, :] = (HCData[i, :] - HC_Predicted[i, :]) / HC_residuals_std_list

    for i in range(TAUData.shape[0]): # 26
        TAU_w[i, :] = (TAUData[i, :] - TAU_Predicted[i, :]) / HC_residuals_std_list_TAU

    for i in range(TDPData.shape[0]): # 30
        TDP_w[i, :] = (TDPData[i, :] - TDP_Predicted[i, :]) / HC_residuals_std_list_TDP

    return HC_w, TAU_w, TDP_w
#### NOT USED ####

def drawthicknessboxplot(HCData, TAUData, TDPData, x_label, y_label, title, outputDir, outputName):
    # HCData, TAUData, TDPData --> Shape N x 400
    # Get Mean for each Regions
    HC_Mean = np.mean(HCData, axis=0)
    TAU_Mean = np.mean(TAUData, axis=0)
    TDP_Mean = np.mean(TDPData, axis=0)

    # Define data
    data = [HC_Mean, TAU_Mean, TDP_Mean]

    # Define figure
    fig, ax = plt.subplots()

    # Draw the boxplots with x_label and y_label
    bplot = ax.boxplot(data, notch=True, labels=x_label)
    ax.set_ylabel(y_label)

    # perform t-test
    t_stat1, pval1 = stats.ttest_ind(HC_Mean, TAU_Mean)
    t_stat2, pval2 = stats.ttest_ind(HC_Mean, TDP_Mean)
    t_stat3, pval3 = stats.ttest_ind(TAU_Mean, TDP_Mean)

    # set title
    ax.set_title(title + f", p(HC vs TAU)={pval1}, p(HC vs TDP)={pval2}, p(TAU vs TDP)={pval3}", fontsize = 5)

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='tif')

def save3DGraph(adjMtx, outputDir, filename):

    # Substitute NaN values to 0
    adjMtx[np.isnan(adjMtx)] = 0

    # Get upper triangle of adjMtx
    adjMtx = np.triu(adjMtx)

    # Get graph of adjMtx
    G = nx.Graph(adjMtx)

    # save graph object to file
    pickle.dump(G, open(outputDir + '/Graph_Objects/' + filename + '.pickle', 'wb'))


def fixedDensity(covMat, N):

    # Find the smallest Nth value. 
    flatArray = covMat.flatten()

    # Sort the array
    flatArray_Sorted = np.sort(flatArray)

    # Get the Nth smallest value (P-value)
    threshVal = flatArray_Sorted[N-1]

    # Return a CovMat, that only keeps the smallest N values. Everything else --> NaN
    fixedDcovMat = np.where(covMat < threshVal, covMat, np.nan) # True / False

    return fixedDcovMat