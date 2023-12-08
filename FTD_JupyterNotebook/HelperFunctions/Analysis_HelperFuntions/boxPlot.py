import matplotlib.pyplot as plt
import scipy
import numpy as np
import os
from scipy import stats
import seaborn as sns
from PIL import Image

def drawthicknessboxplot(HCData, TAUData, TDPData, x_label, y_label, title, outputDir, outputName):
    # HCData, TAUData, TDPData --> Shape N x 400
    # Get Mean for each Regions
    HC_Mean = np.nanmean(HCData, axis=0)
    TAU_Mean = np.nanmean(TAUData, axis=0)
    TDP_Mean = np.nanmean(TDPData, axis=0)

    # Define data
    data = [HC_Mean, TAU_Mean, TDP_Mean]

    # Define figure
    fig, ax = plt.subplots()

    # Draw the boxplots with x_label and y_label
    bplot = ax.boxplot(data, notch=False, labels=x_label)
    ax.set_ylabel(y_label)

    # perform t-test
    t_stat1, pval1 = stats.ttest_ind(HC_Mean, TAU_Mean)
    t_stat2, pval2 = stats.ttest_ind(HC_Mean, TDP_Mean)
    t_stat3, pval3 = stats.ttest_ind(TAU_Mean, TDP_Mean)

    # set title
    ax.set_title(title + f", p(HC vs TAU)={pval1}, p(HC vs TDP)={pval2}, p(TAU vs TDP)={pval3}", fontsize = 5)

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='png')