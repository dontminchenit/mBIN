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
  * Directory Path Setup.ipynb
  * [Data](./FTD_JupyterNotebook/Data)
    * [Pathology_Data](./FTD_JupyterNotebook/Data/Pathology_Data)
    * [ThicknessAtPath_Data](./FTD_JupyterNotebook/Data/ThicknessAtPath_Data)
    * [Thickness_Data](./FTD_JupyterNotebook/Data/Thickness_Data)
  * [HelperFunctions](./FTD_JupyterNotebook/HelperFunctions)
    * [Analysis_HelperFuntions](./FTD_JupyterNotebook/HelperFunctions/Analysis_HelperFuntions)
    * [Load_Dataset_HelperFunctions](./FTD_JupyterNotebook/HelperFunctions/Load_Dataset_HelperFunctions)
  * [Load_Dataset](./FTD_JupyterNotebook/Load_Dataset)
  * [Pathology_Analysis](./FTD_JupyterNotebook/Pathology_Analysis)
    * [Calculated_Data](./FTD_JupyterNotebook/Pathology_Analysis/Calculated_Data)
    * [Figures](./FTD_JupyterNotebook/Pathology_Analysis/Figures)
  * [Thickness_Analysis](./FTD_JupyterNotebook/Thickness_Analysis)
    * [Calculated_Data](./FTD_JupyterNotebook/Thickness_Analysis/Calculated_Data)
    * [Figures](./FTD_JupyterNotebook/Thickness_Analysis/Figures)
  * [Thickness_At_Pathology_Analysis](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis)
    * [Calculated_Data](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/Calculated_Data)
    * [Figures](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/Figures)
  
## Instructions
### Step 1
Open [Directory Path Setup.ipynb] and set directory paths.
### Step 2
Go to [Load_Dataset](./FTD_JupyterNotebook/Load_Dataset) directory and run the notebooks. This will fetch the raw data from your local computer and format the data for analysis. 
Formatted Data would be saved to [Data](./FTD_JupyterNotebook/Data) directory.
### Step 3
Go to [Pathology_Analysis](./FTD_JupyterNotebook/Pathology_Analysis) direcotry to run experiments for Pathology Data.
Run the notebooks in order (numbered [1], [2], ...) since it would save the calculated data into [Calculated_Data](./FTD_JupyterNotebook/Pathology_Analysis/Calculated_Data).
The figures generated from experiments would be saved into [Figures](./FTD_JupyterNotebook/Pathology_Analysis/Figures)
### Step 4
This is similary applied to [Thickness_Analysis](./FTD_JupyterNotebook/Thickness_Analysis) and [Thickness_At_Pathology_Analysis](./FTD_JupyterNotebook/Thickness_At_Pathology_Analysis)
