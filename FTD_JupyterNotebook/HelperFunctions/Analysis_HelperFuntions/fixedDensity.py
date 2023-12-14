import numpy as np

def fixedDensity(covMat, N):
    # Retunr a upper triangle of covMat with smallest N values
    
    # Substitute NaN values in covMat to 0
    covMat[np.isnan(covMat)] = 0

    # Get upper triangle of covMat
    covMat = np.triu(covMat)

    # Get the non-zero part (used for drawing the edge)
    relevant_elements = covMat[covMat != 0]

    # Sort and get the smallest n values
    elements_extract = np.sort(relevant_elements)[:N]

    # Create a boolean mask where True corresponds to the values in the list
    mask = np.isin(covMat, elements_extract)

    # Set everything else to zero (keep only the smallest n values)
    covMat[~mask] = 0
    
    return covMat
    
    
    
#     # Find the smallest Nth value. 
#     flatArray = covMat.flatten()

#     # Sort the array
#     flatArray_Sorted = np.sort(flatArray)

#     # Get the Nth smallest value (P-value)
#     N = 2 * N # Because it is adjacency matrix
#     threshVal = flatArray_Sorted[N-1]

#     # Return a CovMat, that only keeps the smallest N values. Everything else --> NaN
#     fixedDcovMat = np.where(covMat <= threshVal, covMat, np.nan) # True / False

#     return fixedDcovMat
