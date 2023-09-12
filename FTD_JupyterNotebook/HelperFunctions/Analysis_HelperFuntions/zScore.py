import numpy as np

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