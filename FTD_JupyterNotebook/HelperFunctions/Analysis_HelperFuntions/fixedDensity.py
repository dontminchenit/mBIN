import numpy as np

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