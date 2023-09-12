import numpy as np

def pathObsThresh(DataX, obs_thresh):
    """Find the list Region index where there are less than pre-defined observation count 

    Args:
        DataX (ndarray): log %AO of X Pathology Data
        obs_thresh (int): PreDefined observation count. Any regions' observation less than this value would be dicarded
    
    Returns:
        X_missing_index (int list): list of index from the LabelNames that we are not mapping to 3D Atlas due to lacking observations count 
    """
     # Get list containing number of observed subjects per region 
    obs_X = []

    for i in range(DataX.shape[1]): # iterate over the rows of the 2d array / 40
        non_nan_count = np.count_nonzero(~np.isnan(DataX[:, i])) # Number of Non-NaN in this column
        obs_X.append(non_nan_count)

    # TO numpy array
    obs_X = np.array(obs_X)   
    
    X_missing_index = np.argwhere(obs_X < obs_thresh).flatten()

    return X_missing_index