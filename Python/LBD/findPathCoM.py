import numpy as np
import pandas as pd

def findPathCoM(pathLUT, AtlasToPathLUT, AtlasCoM):
    """
    a function to find Center of Mass for Pathology

    Attributes
    ----------
    pathLUT : 
        
    AtlasToPathLUT : 
        
    AtlasCoM : 
       The Center of Mass for Atlas
    
    Returns
    ----------
    pathCoM : [NL x 3 x 2]
        Mean of Center of Mass (?)
    pathToAtlasIndex :

    """
    # Get number of anatomical regions to match to Atlas
    NL = len(AtlasToPathLUT)

    # numpy to save the Center of Mass of the pathology (anatomical region)
    pathCoM = np.zeros((NL, 3, 2))

    # List of matching anatomical region to atlas index
    pathToAtlasIndex = [[None, None] for _ in range(NL)]

    hemi = ['L', 'R']
    
    for l in range(NL):
        for r in range(2):
            # Get the matched Atlas region
            AOName = AtlasToPathLUT.Matched[l]

            # Get indexes of where pathname match Atlas region and {L, R} hemisphere
            # (modified csv file column name Label.Hemi --> Label_Hemi)
            pathLabelIdx = np.where((pathLUT.PathName.str.match(AOName) & pathLUT.Label_Hemi.str.match(hemi[r])))[0]

            if not (pathLabelIdx.size == 0): # If there is a match
                # Get Label_ID400x7 for the matched rows (modified csv file column name Label.ID400x7 --> Label_ID400x7)
                # Need to offset by 1 (index difference between matlab and python)
                atlasIdx = (pathLUT.Label_ID400x7.loc[pathLabelIdx] - 1).to_numpy()
                
                pathCoM[l,:,r] = np.nanmean(AtlasCoM[atlasIdx,:], axis=0)
                pathToAtlasIndex[l][r] = atlasIdx
    
    return pathCoM, pathToAtlasIndex
