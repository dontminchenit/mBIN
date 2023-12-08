import numpy as np
from sklearn.linear_model import LinearRegression

# def generateWScore(AgeSexHC, AgeSexTAU, AgeSexTDP, N, HCData, TAUData, TDPData):
#     """Generate W Scores for HC/TAU/TDP Data

#     Args:
#         AgeSexHC (ndarray): Age and Sex data of HC 
#         AgeSexTAU (ndarray): Age and Sex data of TAU 
#         AgeSexTDP (ndarray): Age and Sex data of TDP 
#         N (int):  Number of Regions we are analyzing in specified thickness Network / N = len(regNames)
#         HCData (ndarray): HC data
#         TAUData (ndarray): TAU data
#         TDPData (ndarray): TDP data

#     Returns:
#         HC_w (float): HC W Score
#         TAU_w (float): TAU W Score
#         TDP_w (float): TDP W Score
#     """

#     HC_model_list = [] # (400,) 
#     HC_residuals_std_list = [] # (400,)  

#     # Linear regression on HC for ALL 400 Regions
#     for k in range(N): 
#         yHC = HCData[:,k] # Thickness values of specified region for HC (54,)
#         assert(yHC.shape == (HCData.shape[0],))
    
#         # Linear Regression
#         regHC = LinearRegression().fit(AgeSexHC, yHC)

#         # Predict 
#         HC_Predict = regHC.predict(AgeSexHC) # Shape (54, ) --> Thickness values for 54 subject in specified Region
#         assert(HC_Predict.shape == (HCData.shape[0],))

#         # Residual
#         HC_Residual = yHC - HC_Predict # Shape (54,)
#         assert(HC_Residual.shape == (HCData.shape[0],))

#         HC_std = np.nanstd(HC_Residual, axis=0, ddof=0)

#         # Save the Linear Regression Model
#         HC_model_list.append(regHC)

#         # Save residual std
#         HC_residuals_std_list.append(HC_std)

#     # Sanity Check
#     assert(len(HC_model_list) == N) # list length of 400
#     assert(len(HC_residuals_std_list) == N) # list length of 400

#     # Predict values of HC, TAU and TDP using these coef / (54, 400), (26, 400), (30, 400)
#     HC_Predicted = np.empty(HCData.shape) # (54, 400) 
#     for h in range(HCData.shape[0]): # 54
#         for g in range (HCData.shape[1]): # 400 or 40
#             HC_Predicted[h, g] = HC_model_list[g].predict([AgeSexHC[h]])

#     TAU_Predicted = np.empty(TAUData.shape) # (26, 400)
#     for h in range(TAUData.shape[0]): # 26
#         for g in range (TAUData.shape[1]): # 400 or 35
#             TAU_Predicted[h, g] = HC_model_list[g].predict([AgeSexTAU[h]])

#     TDP_Predicted = np.empty(TDPData.shape) # (30, 400) 
#     for h in range(TDPData.shape[0]): # 30
#         for g in range (TDPData.shape[1]): # 400 or 35
#             TDP_Predicted[h, g] = HC_model_list[g].predict([AgeSexTDP[h]])

#     # Compute W Score
#     HC_w = np.empty(HCData.shape) # 54 x 400
#     TAU_w = np.empty(TAUData.shape) # 26 x 400
#     TDP_w = np.empty(TDPData.shape) # 30 x 400

#     for i in range(HCData.shape[0]): # 54
#         HC_w[i, :] = (HCData[i, :] - HC_Predicted[i, :]) / HC_residuals_std_list

#     for i in range(TAUData.shape[0]): # 26
#         TAU_w[i, :] = (TAUData[i, :] - TAU_Predicted[i, :]) / HC_residuals_std_list

#     for i in range(TDPData.shape[0]): # 30
#         TDP_w[i, :] = (TDPData[i, :] - TDP_Predicted[i, :]) / HC_residuals_std_list

#     return HC_w, TAU_w, TDP_w




def generateWScore(hcCov, tauCov, tdpCov, N, HCData, TAUData, TDPData):
    """Generate W Scores for HC/TAU/TDP Data

    Args:
        hcCov (ndarray): Covariates of HC data - Age, Sex, and ICV / Row: Subject & Columns: Covariates
        tauCov (ndarray): Covariates of TAU data - Age, Sex, and ICV / Row: Subject & Columns: Covariates
        tdpCov (ndarray): Covariates of TDP data - Age, Sex, and ICV / Row: Subject & Columns: Covariates
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
        
        # Remove NaNs for Linear Regression
        # Identify indices of non-NaN values in yHC
        non_nan_indices = ~np.isnan(yHC)
        # Remove NaN values from yHC
        no_nan_yHC = yHC[non_nan_indices]

        # Remove corresponding rows from hcCov
        no_nan_hcCov = hcCov[non_nan_indices]
        
        # Linear Regression
        regHC = LinearRegression().fit(no_nan_hcCov, no_nan_yHC)

        # Predict 
        HC_Predict = regHC.predict(hcCov) # Shape (54, ) --> Thickness values for 54 subject in specified Region
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
            HC_Predicted[h, g] = HC_model_list[g].predict([hcCov[h]])

    TAU_Predicted = np.empty(TAUData.shape) # (26, 400)
    for h in range(TAUData.shape[0]): # 26
        for g in range (TAUData.shape[1]): # 400 or 35
            TAU_Predicted[h, g] = HC_model_list[g].predict([tauCov[h]])

    TDP_Predicted = np.empty(TDPData.shape) # (30, 400) 
    for h in range(TDPData.shape[0]): # 30
        for g in range (TDPData.shape[1]): # 400 or 35
            TDP_Predicted[h, g] = HC_model_list[g].predict([tdpCov[h]])

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