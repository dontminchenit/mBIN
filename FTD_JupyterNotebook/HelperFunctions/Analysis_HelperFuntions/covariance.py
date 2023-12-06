import numpy as np
import numpy.ma as ma
from scipy.stats import norm

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

def covCal(x, y, cov_thresh):
    """ Function that calculates Covariance Matrix between x and y
    
    Args:
        x (ndarray): 1st dataset
        y (ndarray): 2nd dataset
        cov_thresh (float): covariance matrix threshold (just to make suere there is no noise)
    
    Returns:
        covMat (ndarray) cov mat we calculated for x and y
    """
    
    # Sanity Check - Type Check
    assert(isinstance(x,np.ndarray))
    assert(isinstance(y,np.ndarray))
    
    # Sanity Check - if dataset have same number of regions (columns)
    assert(x.shape[1] == y.shape[1])
    
    # Define N (Number of regions, columns)
    N = x.shape[1]
    
    # Define Empty Cov Mat
    covMat = np.full((N,N), np.nan)
    
    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of different regions for x and y both not NaN
                n = np.sum(~np.isnan(x[:,i]) & ~np.isnan(y[:,j]))
                
                if n > 3:
                    covMat[i, j] = ma.corrcoef(ma.masked_invalid(x[:,i]), ma.masked_invalid(y[:,j]))[0, 1]
                    if covMat[i, j] < cov_thresh:
                        covMat[i, j] = np.nan
    
    return covMat

def covCalSigXY(x, y, covMatX, covMatY, pthresh, cov_thresh):
    """ Function that calculates Covariance Matrix X>Y and Y>X
    
    Args:
        x (ndarray): 1st dataset
        y (ndarray): 2nd dataset
        covMatX (ndarray): covMat of x vs x
        covMatY (ndarray): covMat of y vs y
        pthresh (float): p-value threshold
        cov_thresh (float): covariance matrix threshold (just to make suere there is no noise)
    
    Returns:
        covMatXgtY (ndarray) cov mat we calculated for x and y, x greater than y.  X > Y
        covMatYgtX (ndarray) cov mat we calculated for x and y, y greater than x.  Y > X
    """
    
    # Sanity Check - Type Check
    assert(isinstance(x,np.ndarray))
    assert(isinstance(y,np.ndarray))
    
    # Sanity Check - if dataset have same number of regions (columns)
    assert(x.shape[1] == y.shape[1])
    
    # Define N (Number of regions, columns)
    N = x.shape[1]
    
    # Define Empty Cov Mat
    covMatXgtY = np.full((N,N), np.nan)
    covMatYgtX = np.full((N,N), np.nan)
    
    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of different regions where x and y both are not NaN
                nX = np.sum(~np.isnan(x[:,i]) & ~np.isnan(x[:,j]))
                nY = np.sum(~np.isnan(y[:,i]) & ~np.isnan(y[:,j]))
                
                if nX > 3:
                    if covMatX[i, j] < cov_thresh:
                        nX = 0
                if nY > 3:
                    if covMatY[i, j] < cov_thresh:
                        nY = 0
                
                if nX > 3 and nY > 3:
                    covMatXgtY[i,j] = test_corr_sig(covMatY[i,j],covMatX[i,j],nY,nX)[0] < pthresh 
                    
                    covMatYgtX[i,j] = test_corr_sig(covMatX[i,j],covMatY[i,j],nX,nY)[0] < pthresh

    return covMatXgtY, covMatYgtX
                
def covCalSigXYRaw(x, y, covMatX, covMatY, cov_thresh):
    """ Function that calculates Covariance Matrix X>Y and Y>X
    
    Args:
        x (ndarray): 1st dataset
        y (ndarray): 2nd dataset
        covMatX (ndarray): covMat of x vs x
        covMatY (ndarray): covMat of y vs y
        pthresh (float): p-value threshold
        cov_thresh (float): covariance matrix threshold (just to make suere there is no noise)
    
    Returns:
        covMatXgtY (ndarray) cov mat we calculated for x and y, x greater than y.  X > Y
        covMatYgtX (ndarray) cov mat we calculated for x and y, y greater than x.  Y > X
    """
    
    # Sanity Check - Type Check
    assert(isinstance(x,np.ndarray))
    assert(isinstance(y,np.ndarray))
    
    # Sanity Check - if dataset have same number of regions (columns)
    assert(x.shape[1] == y.shape[1])
    
    # Define N (Number of regions, columns)
    N = x.shape[1]
    
    # Define Empty Cov Mat
    covMatXgtYRaw = np.full((N,N), np.nan)
    covMatYgtXRaw = np.full((N,N), np.nan)
    
    for i in range(N):
        for j in range(N):
            if i!=j:
                # Sum of different regions where x and y both are not NaN
                nX = np.sum(~np.isnan(x[:,i]) & ~np.isnan(x[:,j]))
                nY = np.sum(~np.isnan(y[:,i]) & ~np.isnan(y[:,j]))
                
                if nX > 3:
                    if covMatX[i, j] < cov_thresh:
                        nX = 0
                if nY > 3:
                    if covMatY[i, j] < cov_thresh:
                        nY = 0
                
                if nX > 3 and nY > 3:
                    covMatXgtYRaw[i,j] = test_corr_sig(covMatY[i,j],covMatX[i,j],nY,nX)[0]
                    
                    covMatYgtXRaw[i,j] = test_corr_sig(covMatX[i,j],covMatY[i,j],nX,nY)[0]

    return covMatXgtYRaw, covMatYgtXRaw                
      
def calCovMatDict(dataHC, dataTAU, dataTDP, pthresh, cov_thresh):
    cov_dict = {}
    
    covMatHC = covCal(dataHC, dataHC, cov_thresh)
    covMatTAU = covCal(dataTAU, dataTAU, cov_thresh)
    covMatTDP = covCal(dataTDP, dataTDP, cov_thresh)
    
    cmpCovTAU_gt_TDP, cmpCovTDP_gt_TAU = covCalSigXY(dataTAU, dataTDP, covMatTAU, covMatTDP, pthresh, cov_thresh)
    cmpCovTAU_gt_HC, cmpCovTAU_lt_HC = covCalSigXY(dataTAU, dataHC, covMatTAU, covMatHC, pthresh, cov_thresh)
    cmpCovTDP_gt_HC, cmpCovTDP_lt_HC = covCalSigXY(dataTDP, dataHC, covMatTDP, covMatHC, pthresh, cov_thresh)
    
    cmpCovTAU_gt_TDP_raw, cmpCovTDP_gt_TAU_raw = covCalSigXYRaw(dataTAU, dataTDP, covMatTAU, covMatTDP, cov_thresh)
    cmpCovTAU_gt_HC_raw, cmpCovTAU_lt_HC_raw = covCalSigXYRaw(dataTAU, dataHC, covMatTAU, covMatHC, cov_thresh)
    cmpCovTDP_gt_HC_raw, cmpCovTDP_lt_HC_raw = covCalSigXYRaw(dataTDP, dataHC, covMatTDP, covMatHC, cov_thresh)
    
    cov_dict["HC"] = covMatHC
    cov_dict["TAU"] = covMatTAU
    cov_dict["TDP"] = covMatTDP
    
    cov_dict["TAU_gt_TDP"] = cmpCovTAU_gt_TDP
    cov_dict["TDP_gt_TAU"] = cmpCovTDP_gt_TAU
    cov_dict["TAU_gt_HC"] = cmpCovTAU_gt_HC
    cov_dict["TAU_lt_HC"] = cmpCovTAU_lt_HC
    cov_dict["TDP_gt_HC"] = cmpCovTDP_gt_HC
    cov_dict["TDP_lt_HC"] = cmpCovTDP_lt_HC
    
    cov_dict["TAU_gt_TDP_raw"] = cmpCovTAU_gt_TDP_raw
    cov_dict["TDP_gt_TAU_raw"] = cmpCovTDP_gt_TAU_raw
    cov_dict["TAU_gt_HC_raw"] = cmpCovTAU_gt_HC_raw
    cov_dict["TAU_lt_HC_raw"] = cmpCovTAU_lt_HC_raw
    cov_dict["TDP_gt_HC_raw"] = cmpCovTDP_gt_HC_raw
    cov_dict["TDP_lt_HC_raw"] = cmpCovTDP_lt_HC_raw
    
    return cov_dict
                
                
                
                
                
                
                
                
                