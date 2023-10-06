import matplotlib.pyplot as plt
import scipy
import numpy as np
import os
from scipy import stats
import seaborn as sns
from PIL import Image

# Correlation between Nodal Strength vs Actual Value
def nonZeroDegCorr(DataX, covMatX, ymin, ymax, title, x_label, y_label, outputDir, outputName, linear_regression = False):
    
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
   
    # set xaxes range
    plt.ylim(ymin, ymax)
    
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

def nonZeroDegCorrCloseFar(DataX, covMatX, close_connection_list, 
                           ymin, ymax, title, x_label, y_label, 
                           outputDir, outputName, linear_regression = False, savefig=True):
    
    # Copy the Covariance Matrix and set negative values as zero
    covMatXnz = covMatX.copy()
    covMatXnz[covMatXnz < 0] = 0
    
    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degX = np.sum(covMatXnz, axis=0, where=~np.isnan(covMatXnz)) 
    
    # Divide the Covariance Matrix to Close and Far
    deg_close_array = []
    for i in range(covMatXnz.shape[1]): # For each column in array
        # Get the column
        col = covMatXnz[:, i]

        # Get the Close indices for that col (region)
        close_ind = close_connection_list[i]

        # Get the Degree for Close
        deg_close = np.nansum(col[close_ind])
        deg_close_array.append(deg_close)

    # Convert Deg Close to numpy
    deg_close = np.array(deg_close_array)
    
    # Get Deg Far
    deg_far = degX - deg_close_array
    
    # Normalize the Degree Sum (Due to difference in number of nodes that is close vs far)
    # Get Number of Close/Far Nodes
    closeNum = [len(sublist) for sublist in close_connection_list]
    farNum = (np.zeros((len(close_connection_list), ), dtype=int) + len(close_connection_list)) - closeNum
    
    deg_close = deg_close/closeNum
    deg_far = deg_far/farNum
    
    # Define figure
    fig = plt.figure()
    
    # Figure size
    plt.figure(figsize=(3.5,8))

    # Draw Scatter Plot 
    # Close - empty circle with green edges
    plt.scatter(deg_close, np.nanmean(DataX, axis=0), facecolors='none', edgecolors='orange')
    # Far - empty circle with blue edges
    plt.scatter(deg_far, np.nanmean(DataX, axis=0), facecolors='none', edgecolors='blue')
    # set yaxis range
    plt.ylim(ymin, ymax)
    
    # Get r and p-value
    r_close, p_close = scipy.stats.pearsonr(deg_close, np.nanmean(DataX, axis=0))
    r_far, p_far = scipy.stats.pearsonr(deg_far, np.nanmean(DataX, axis=0))
 
    # Set title
    plt.title(title + f",Close: r={r_close:.6f}, p={p_close:.6f}/Far: r={r_far:.6f}, p={p_far:.6f}", fontsize=6)
   
    # Draw Linear Regression Line (is set to True)
    if linear_regression:
        # Obtain m (slope) and b(intercept) of linear regression line
        m_close, b_close = np.polyfit(deg_close, np.nanmean(DataX, axis=0), 1)
        m_far, b_far = np.polyfit(deg_far, np.nanmean(DataX, axis=0), 1)
        
        #add linear regression line to scatterplot 
        plt.plot(deg_close, m_close*deg_close+b_close, color="red")
        plt.plot(deg_far, m_far*deg_far+b_far, color="green")

    # Set X and Y Labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    if savefig:
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

    
    
    
    
#--------------------
def nonZeroDegCorrCloseFarInv(DataX, covMatX, close_connection_list, 
                           ymin, ymax, title, x_label, y_label, 
                           outputDir, outputName, linear_regression = False):
    
    # Copy the Covariance Matrix and set negative values as zero
    covMatXnz = covMatX.copy()
    covMatXnz[covMatXnz < 0] = 0
    
    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degX = np.sum(covMatXnz, axis=0, where=~np.isnan(covMatXnz)) 
    
    # Divide the Covariance Matrix to Close and Far
    deg_close_array = []
    for i in range(covMatXnz.shape[1]): # For each column in array
        # Get the column
        col = covMatXnz[:, i]

        # Get the Close indices for that col (region)
        close_ind = close_connection_list[i]

        # Get the Degree for Close
        deg_close = np.nansum(col[close_ind])
        deg_close_array.append(deg_close)

    # Convert Deg Close to numpy
    deg_close = np.array(deg_close_array)
    
    # Get Deg Far
    deg_far = degX - deg_close_array
    
    # Define figure
    fig = plt.figure()
    
    # Figure size
    #plt.figure(figsize=(3.5,8))

    # Draw Scatter Plot 
    # Close - empty circle with green edges
    plt.scatter(np.nanmean(DataX, axis=0), deg_close, facecolors='none', edgecolors='orange')
    # Far - empty circle with blue edges
    plt.scatter(np.nanmean(DataX, axis=0), deg_far, facecolors='none', edgecolors='blue')
    # set yaxis range
    plt.xlim(ymin, ymax)
    
    
    
    # Get r and p-value
    r_close, p_close = scipy.stats.pearsonr(np.nanmean(DataX, axis=0), deg_close)
    r_far, p_far = scipy.stats.pearsonr(np.nanmean(DataX, axis=0), deg_far)
 
    # Set title
    plt.title(title + f",Close: r={r_close:.6f}, p={p_close:.6f}/Far: r={r_far:.6f}, p={p_far:.6f}", fontsize=6)
   
    # Draw Linear Regression Line (is set to True)
    if linear_regression:
        # Obtain m (slope) and b(intercept) of linear regression line
        m_close, b_close = np.polyfit(np.nanmean(DataX, axis=0), deg_close, 1)
        m_far, b_far = np.polyfit(np.nanmean(DataX, axis=0), deg_far, 1)
        
        #add linear regression line to scatterplot 
        plt.plot(np.nanmean(DataX, axis=0), m_close*np.nanmean(DataX, axis=0)+b_close, color="red")
        plt.plot(np.nanmean(DataX, axis=0), m_far*np.nanmean(DataX, axis=0)+b_far, color="green")

    # Set X and Y Labels
    plt.ylabel(x_label)
    plt.xlabel(y_label)

    # Save Figure
    plt.savefig(os.path.join(outputDir, outputName) + '.png', dpi=400)
    
    # Show Figure
    plt.show()