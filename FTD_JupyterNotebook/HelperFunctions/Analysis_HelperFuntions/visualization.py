import matplotlib.pyplot as plt
import scipy
import numpy as np
import os
from scipy import stats
import seaborn as sns
from PIL import Image
from chordPlot import condition

def drawCovMatrix(covMat, x_labels, y_labels, title, outputDir, annot_fontsize = 8, tick_fontsize = 6, annot_bool = True, save=False):
    # Dimension(one side) of the covariance matrix
    N = covMat.shape[0]

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(covMat, dtype=bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots()

    # colormap - 'cold' 'hot'
    # cmap = sns.color_palette("coolwarm", as_cmap=True)
    cmap = 'Reds'
    
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(covMat, mask=mask, cmap=cmap, vmax=1, vmin=0, center=0.5,
                square=True, linewidths=.5, xticklabels=x_labels, yticklabels=y_labels, cbar_kws={"shrink": .5}, 
                annot=annot_bool, annot_kws={"size": annot_fontsize}, fmt='.2f')
    
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)

    plt.tight_layout()

    # Set figure title
    plt.title(title)
    
    # Save
    if save:
        plt.savefig(os.path.join(outputDir, title + '.png'), dpi=400)
    
    # Show the figure
    plt.show()
    
def drawCovMatrixSubpath(covMat, x_labels, y_labels, title, outputDir, subpath, annot_fontsize=8, tick_fontsize=6, annot_bool=True, save=False):
    N = covMat.shape[0]
    mask = np.triu(np.ones_like(covMat, dtype=bool))
    
    # New mask based on the condition
    condition_mask = np.zeros_like(covMat, dtype=bool)
    
    # Apply condition to filter the connections
    for i, source in enumerate(x_labels):
        for j, target in enumerate(y_labels):
            if condition(source, target, subpath):
                condition_mask[i, j] = True

    
    # Combine the condition mask with the upper triangle mask
    mask = mask | ~condition_mask.T

    f, ax = plt.subplots()
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    sns.heatmap(covMat, mask=mask, cmap=cmap, vmax=1, vmin=-1, center=0,
                square=True, linewidths=.5, xticklabels=x_labels, yticklabels=y_labels, cbar_kws={"shrink": .5}, 
                annot=annot_bool, annot_kws={"size": annot_fontsize}, fmt='.2f')
    
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.title(title)
    
    if save:
        plt.savefig(os.path.join(outputDir, title + '.png'), dpi=400)
    
    plt.show()

def drawCovMatrixStrat(covMat, x_labels, y_labels, title, outputDir, minval, maxval, annot_fontsize=8, tick_fontsize=6, annot_bool=True, save=False):
    N = covMat.shape[0]
    mask = np.triu(np.ones_like(covMat, dtype=bool))
    
    # New mask based on the condition
    condition_mask = np.zeros_like(covMat, dtype=bool)
    
    # Apply condition to filter the connections
    for i, source in enumerate(x_labels):
        for j, target in enumerate(y_labels):
            if (covMat[i, j] >= minval) & (covMat[i, j] < maxval):
                condition_mask[i, j] = True

    
    # Combine the condition mask with the upper triangle mask
    mask = mask | ~condition_mask.T

    f, ax = plt.subplots()
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    sns.heatmap(covMat, mask=mask, cmap=cmap, vmax=1, vmin=-1, center=0,
                square=True, linewidths=.5, xticklabels=x_labels, yticklabels=y_labels, cbar_kws={"shrink": .5}, 
                annot=annot_bool, annot_kws={"size": annot_fontsize}, fmt='.2f')
    
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.title(title)
    
    if save:
        plt.savefig(os.path.join(outputDir, title + '.png'), dpi=400)
    
    plt.show()