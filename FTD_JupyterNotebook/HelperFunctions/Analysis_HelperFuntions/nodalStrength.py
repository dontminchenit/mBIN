import matplotlib.pyplot as plt
import scipy
import numpy as np
import os
from scipy import stats
import seaborn as sns
from PIL import Image

# Correlation between Nodal Strength vs Actual Value
def nonZeroDegCorr(DataX, covMatX, title, x_label, y_label, outputDir, outputName, linear_regression = False):
    
    # Copy the Covariance Matrix and set negative values as zero
    covMatXnz = covMatX.copy()
    covMatXnz[covMatXnz < 0] = 0

    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degX = np.sum(covMatXnz, axis=0, where=~np.isnan(covMatXnz)) 

    # Define figure
    fig = plt.figure()
    
    # Figure size
    plt.figure(figsize=(3.5,8))

    # Draw Scatter Plot / empty circle with blue edges
    plt.scatter(degX, np.nanmean(DataX, axis=0), facecolors='none', edgecolors='b')
   
    # Get r and p-value
    r, p = scipy.stats.pearsonr(degX, np.nanmean(DataX, axis=0))
 
    # Set title
    plt.title(title + f", r={r:.6f}, p={p:.6f}", fontsize=10)
   
    # Draw Linear Regression Line (is set to True)
    if linear_regression:
        # Obtain m (slope) and b(intercept) of linear regression line
        m, b = np.polyfit(degX, np.nanmean(DataX, axis=0), 1)
        
        #add linear regression line to scatterplot 
        plt.plot(degX, m*degX+b, color="red")

    # Set X and Y Labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Save Figure
    plt.savefig(os.path.join(outputDir, outputName) + '.png', dpi=400)
    
    # Show Figure
    plt.show()
    
# Compare Distribution(Boxplot) of Nodal Strength Between HC, TAU, and TDP
def nodalStrengthDist(covMatA, covMatB, covMatC, x_label, y_label, title, outputDir, outputName):
    
    # Copy the Covariance Matrix and set negative values as zero.
    covMatAnz = covMatA.copy()
    covMatAnz[covMatAnz < 0] = 0
    
    covMatBnz = covMatB.copy()
    covMatBnz[covMatBnz < 0] = 0
    
    covMatCnz = covMatC.copy()
    covMatCnz[covMatCnz < 0] = 0

    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degA = np.sum(covMatAnz, axis=0, where=~np.isnan(covMatAnz)) # Shape: (91, )
    degB = np.sum(covMatBnz, axis=0, where=~np.isnan(covMatBnz)) # Shape: (91, )
    degC = np.sum(covMatCnz, axis=0, where=~np.isnan(covMatCnz)) # Shape: (91, )

    # Define data
    data = [degA, degB, degC]

    # Define figure
    fig, ax = plt.subplots()

    # Draw the boxplots with x_label and y_label
    bplot = ax.boxplot(data, notch=True, labels=x_label)
    ax.set_ylabel(y_label)

    # perform t-test
    t_stat, pval = stats.ttest_ind(degB, degC)

    # set title
    ax.set_title(title + f", p_value between TAU vs TDP p={pval}")

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='png')