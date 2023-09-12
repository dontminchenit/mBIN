# Multimodal Brain Integration Networks (mBIN)
Multimodal Brain Integration Networks (mBIN) is toolbox for high dimensional analysis of brain networks across multiple imaging modalities. Our goal is to enable high level visualization and interpretation of  structual and functional brain data across micro to macro scales. 

This repository includes analysis software used in the following publications:

Chen, M., Ohm, D.T., Phillips, J.S., McMillan, C.T., Capp, N., Peterson, C., Xie, E., Wolk, D.A., Trojanowski, J.Q., Lee, E.B., Gee, J., Grossman, M. and Irwin, D.J., 2022. Divergent histopathological networks of frontotemporal degeneration proteinopathy subytpes. Journal of Neuroscience, 42(18), pp.3868-3877.

Chen, M., Burke, S., Olm, C.A., Irwin, D.J., Massimo, L., Lee, E.B., Trojanowski, J.Q, Gee, J.C, Grossman, M., 2023, Antemortem Network Analysis of Spreading Pathology in Autopsy-Confirmed Frontotemporal Degeneration. Brain Communications. Accepted.

FAQ:
Version Compatibility: 
  You need "matlab statistics and machine learning toolbox" to run the code;
  You need "Matlab 2022" to run the function "clim" in "Matlab\NetworkAnalysisGeneral\Code\SupportFunctions\plotgiiSurf.m"

# Instructions for Running Jupyter Notebook
## Directory Structure
* [FTD_JupyterNotebook](./FTD_JupyterNotebook)
  * [Data/](./FTD_JupyterNotebook/Data)
    * [Pathology_Data/](./FTD_JupyterNotebook/Data/Pathology_Data)
    * [ThicknessAtPath_Data/](./FTD_JupyterNotebook/Data/ThicknessAtPath_Data)
    * [Thickness_Data/](./FTD_JupyterNotebook/Data/Thickness_Data)
  * [HelperFunctions/](./FTD_JupyterNotebook/HelperFunctions)
    * [Analysis_HelperFuntions/](./FTD_JupyterNotebook/HelperFunctions/Analysis_HelperFuntions)
    * [Load_Dataset_HelperFunctions/](./FTD_JupyterNotebook/HelperFunctions/Load_Dataset_HelperFunctions)
  * [Load_Dataset/](./FTD_JupyterNotebook/Load_Dataset)
    * [Load_Data-Pathology.ipynb](./FTD_JupyterNotebook/Load_Dataset/Load_Data-Pathology.ipynb)
    * [Load_Data-Thickness.ipynb](./FTD_JupyterNotebook/Load_Dataset/Load_Data-Thickness.ipynb)
    * [Load_Data-ThicknessAtPath.ipynb](./FTD_JupyterNotebook/Load_Dataset/Load_Data-ThicknessAtPath.ipynb)
  * [Pathology_Analysis/](./FTD_JupyterNotebook/Pathology_Analysis)
    * [.ipynb_checkpoints/](./FTD_JupyterNotebook/Pathology_Analysis/.ipynb_checkpoints)
    * [Calculated_Data/](./FTD_JupyterNotebook/Pathology_Analysis/Calculated_Data)
    * [Figures/](./FTD_JupyterNotebook/Pathology_Analysis/Figures)
    * [[1] Compute Pathology Covariance Matrix [GM].ipynb](./FTD_JupyterNotebook/Pathology_Analysis/[1] Compute Pathology Covariance Matrix [GM].ipynb)
    * [[2] Drop Regions with Few Observations (Outliers) [GM].ipynb](./FTD_JupyterNotebook/Pathology_Analysis/[2] Drop Regions with Few Observations (Outliers) [GM].ipynb)
    * [[3] Covariance Matrices Visualization.ipynb](./FTD_JupyterNotebook/Pathology_Analysis/[3] Covariance Matrices Visualization.ipynb)
    * [[4] Nodal Strength vs Log %AO.ipynb](./FTD_JupyterNotebook/Pathology_Analysis/[4] Nodal Strength vs Log %AO.ipynb)
    * [[5] 3D Atlas Mapping.ipynb](./FTD_JupyterNotebook/Pathology_Analysis/[5] 3D Atlas Mapping.ipynb)
  * [Thickness_Analysis/](./FTD_JupyterNotebook/Thickness_Analysis)
    * [.ipynb_checkpoints/](./FTD_JupyterNotebook/Thickness_Analysis/.ipynb_checkpoints)
    * [Calculated_Data/](./FTD_JupyterNotebook/Thickness_Analysis/Calculated_Data)
    * [Figures/](./FTD_JupyterNotebook/Thickness_Analysis/Figures)
    * [[1] Network Selection.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[1] Network Selection.ipynb)
    * [[2] Z, W-Score Calculation.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[2] Z, W-Score Calculation.ipynb)
    * [[3] Compute Thickness Covariance Matrix.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[3] Compute Thickness Covariance Matrix.ipynb)
    * [[4] Covariance Matrix Visualization.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[4] Covariance Matrix Visualization.ipynb)
    * [[5] Thickness Value (Original, Z Score, W Score) Distribution.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[5] Thickness Value (Original, Z Score, W Score) Distribution.ipynb)
    * [[6] Normalized Volume Value (Original, Z Score, W Score) Distribution.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[6] Normalized Volume Value (Original, Z Score, W Score) Distribution.ipynb)
    * [[7] Nodal Strength vs Thickness.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[7] Nodal Strength vs Thickness.ipynb)
    * [[8] Nodal Strength Distribution.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[8] Nodal Strength Distribution.ipynb)
    * [[9] 3D Atlas Mapping.ipynb](./FTD_JupyterNotebook/Thickness_Analysis/[9] 3D Atlas Mapping.ipynb)
  * [Thickness_At_Pathology_Analysis/](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis)
    * [.ipynb_checkpoints/](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/.ipynb_checkpoints)
    * [Calculated_Data/](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/Calculated_Data)
    * [Figures/](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/Figures)
    * [[1] Z, W-Score Calculation.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[1] Z, W-Score Calculation.ipynb)
    * [[2] Compute Thickness At Path Covariance Matrix.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[2] Compute Thickness At Path Covariance Matrix.ipynb)
    * [[3] Drop Regions with Few Observations (Outliers).ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[3] Drop Regions with Few Observations (Outliers).ipynb)
    * [[4] Thickness Value At Pathology (Original, Z Score, W Score) Distribution.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[4] Thickness Value At Pathology (Original, Z Score, W Score) Distribution.ipynb)
    * [[5] Normalized Volume At Pathology (Original, Z Score, W Score) Distribution.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[5] Normalized Volume At Pathology (Original, Z Score, W Score) Distribution.ipynb)
    * [[6] Nodal Strength vs Thickness At Pathology.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[6] Nodal Strength vs Thickness At Pathology.ipynb)
    * [[7] Nodal Strength Distribution.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[7] Nodal Strength Distribution.ipynb)
    * [[8] 3D Atlas Mapping.ipynb](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/[8] 3D Atlas Mapping.ipynb)
  * [Directory Path Setup.ipynb](./FTD_JupyterNotebook/Directory Path Setup.ipynb)
