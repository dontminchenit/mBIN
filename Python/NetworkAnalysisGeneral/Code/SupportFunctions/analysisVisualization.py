import matplotlib.pyplot as plt
import scipy
import numpy as np
import os
from scipy import stats
import seaborn as sns
from PIL import Image

def drawScatterplot(x, y, x_label, y_label, title, outputDir, outputName, linear_regression = False):
    # create scatterplot / empty circle with blue edges
    plt.scatter(x, y, facecolors='none', edgecolors='b')
    
    # set x and y axis label
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Pearson's r and p-value
    r, p = scipy.stats.pearsonr(x, y)    

    # Set figure title / including r and p-value
    plt.title(title + f", r={r}, p={p}")

    if linear_regression:
        # Obtain m (slope) and b(intercept) of linear regression line
        m, b = np.polyfit(x, y, 1)

        #add linear regression line to scatterplot 
        plt.plot(x, m*x+b)
    
    # Save the figure
    plt.savefig(os.path.join(outputDir, outputName), dpi=400)

def nonZeroDegCmp(covMatCtrl, covMatyesAD, covMatnoAD, x_label, y_label, title, outputDir, outputName):
    
    # Copy the Covariance Matrix and set negative values as zero.
    covMatCtrlnz = covMatCtrl.copy()
    covMatCtrlnz[covMatCtrlnz < 0] = 0
    
    covMatyesADnz = covMatyesAD.copy()
    covMatyesADnz[covMatyesADnz < 0] = 0

    covMatnoADnz = covMatnoAD.copy()
    covMatnoADnz[covMatnoADnz < 0] = 0

    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degYes = np.sum(covMatyesADnz, axis=0, where=~np.isnan(covMatyesADnz)) # Shape: (91, )
    degNo = np.sum(covMatnoADnz, axis=0, where=~np.isnan(covMatnoADnz)) # Shape: (91, )
    degCtrl = np.sum(covMatCtrlnz, axis=0, where=~np.isnan(covMatCtrlnz)) # Shape: (91, )

    # Define data
    data = [degCtrl, degYes, degNo]

    # Define figure
    fig, ax = plt.subplots()

    # Draw the boxplots with x_label and y_label
    bplot = ax.boxplot(data, notch=True, labels=x_label)
    ax.set_ylabel(y_label)

    # perform t-test
    t_stat, pval = stats.ttest_ind(degYes, degNo)

    # set title
    ax.set_title(title + f", p={pval}")

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='tif')

def nonZeroDegCorr(ctrl, thickyesAD, thicknoAD, covMatCtrl, covMatyesAD, covMatnoAD, subplot1_title, subplot2_title, subplot3_title, x_label, y_label, outputDir, outputName, linear_regression = False):
    
    # Copy the Covariance Matrix and set negative values as zero.
    covMatCtrlnz = covMatCtrl.copy()
    covMatCtrlnz[covMatCtrlnz < 0] = 0

    covMatyesADnz = covMatyesAD.copy()
    covMatyesADnz[covMatyesADnz < 0] = 0

    covMatnoADnz = covMatnoAD.copy()
    covMatnoADnz[covMatnoADnz < 0] = 0

    # Get sum of covariance values for all regions respective to each region
    # Similar to computing the degree of nodes in a Network
    degCtrl = np.sum(covMatCtrlnz, axis=0, where=~np.isnan(covMatCtrlnz)) 
    degYes = np.sum(covMatyesADnz, axis=0, where=~np.isnan(covMatyesADnz)) # Shape: (91, )
    degNo = np.sum(covMatnoADnz, axis=0, where=~np.isnan(covMatnoADnz)) # Shape: (91, )

    # Define figure
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)


    # Draw Scatter Plot / empty circle with blue edges
    ax1.scatter(degCtrl, np.mean(ctrl, axis=0), facecolors='none', edgecolors='b')
    ax2.scatter(degYes, np.mean(thickyesAD, axis=0), facecolors='none', edgecolors='b') # Degree [1 x 91] vs mean of thickness value of Patients [1 x 91]
    ax3.scatter(degNo, np.mean(thicknoAD, axis=0), facecolors='none', edgecolors='b') # Degree [1 x 91] vs mean of thickness value of Patients [1 x 91]

    # Get r and p-value
    r1, p1 = scipy.stats.pearsonr(degCtrl, np.mean(ctrl, axis=0))
    r2, p2 = scipy.stats.pearsonr(degYes, np.mean(thickyesAD, axis=0))
    r3, p3 = scipy.stats.pearsonr(degNo, np.mean(thicknoAD, axis=0))

    # Set title for each subplot
    ax1.set_title(subplot1_title + f", r={r1:.6f}, p={p1:.6f}", fontsize=5)
    ax2.set_title(subplot2_title + f", r={r2:.6f}, p={p2:.6f}", fontsize=5)
    ax3.set_title(subplot3_title + f", r={r3:.6f}, p={p3:.6f}", fontsize=5)

    # Draw Linear Regression Line (is set to True)
    if linear_regression:
        # Obtain m (slope) and b(intercept) of linear regression line
        m1, b1 = np.polyfit(degCtrl, np.mean(ctrl, axis=0), 1)
        m2, b2 = np.polyfit(degYes, np.mean(thickyesAD, axis=0), 1)
        m3, b3 = np.polyfit(degNo, np.mean(thicknoAD, axis=0), 1)
        

        #add linear regression line to scatterplot 
        ax1.plot(degCtrl, m1*degCtrl+b1, color="red")
        ax2.plot(degYes, m2*degYes+b2, color="red")
        ax3.plot(degNo, m3*degNo+b3, color="red")

    # Set X and Y Labels
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    # ax2.set_ylabel(y_label)
    ax3.set_xlabel(x_label)
    # ax3.set_ylabel(y_label)

    # Set tick fontsize
    ax1.tick_params(axis='both', which='major', labelsize=5)
    ax2.tick_params(axis='both', which='major', labelsize=5)
    ax3.tick_params(axis='both', which='major', labelsize=5)

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='tif')


def drawCovMatrix(covMat, x_labels, y_labels, title, outputDir, outputName, annot_fontsize = 8, tick_fontsize = 6, annot_bool = True):
    # Dimension(one side) of the covariance matrix
    N = covMat.shape[0]

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(covMat, dtype=bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots()

    # colormap - 'cold' 'hot'
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(covMat, mask=mask, cmap=cmap, vmax=1, vmin=0, center=0.5,
                square=True, linewidths=.5, xticklabels=x_labels, yticklabels=y_labels, cbar_kws={"shrink": .5}, annot=annot_bool, annot_kws={"size": annot_fontsize}, fmt='.2f')
    
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)

    # Set figure title
    plt.title(title)

    # Save the figure
    plt.savefig(os.path.join(outputDir, outputName), dpi=400)

def get_concat_h(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def cropImg4(img_list):
    n = len(img_list)

    assert(n == 4)

    # y_top1 = 630
    # y_bot1 = 1720
    # y_top2 = 2800
    # y_bot2 = 4175

    y_top1 = 550
    y_bot1 = 1850
    y_top2 = 2770
    y_bot2 = 4220

    img1_cropped1 = img_list[0].crop((530, y_top1, 1900, y_bot1))
    img1_cropped2 = img_list[0].crop((530, y_top2, 1900, y_bot2))

    img2_cropped1 = img_list[1].crop((530, y_top1, 1900, y_bot1))
    img2_cropped2 = img_list[1].crop((530, y_top2, 1900, y_bot2))

    img3_cropped1 = img_list[2].crop((530, y_top1, 1900, y_bot1))
    img3_cropped2 = img_list[2].crop((530, y_top2, 1900, y_bot2))

    img4_cropped1 = img_list[3].crop((530, y_top1, 1900, y_bot1))
    img4_cropped2 = img_list[3].crop((530, y_top2, 1900, y_bot2))


    img1_comb = get_concat_v(img1_cropped1, img1_cropped2)
    img2_comb = get_concat_v(img2_cropped1, img2_cropped2)
    img3_comb = get_concat_v(img3_cropped1, img3_cropped2)
    img4_comb = get_concat_v(img4_cropped1, img4_cropped2)

    img12 = get_concat_h(img1_comb, img2_comb)
    img34 = get_concat_h(img3_comb, img4_comb)

    final_img = get_concat_h(img12, img34)
    return final_img

def cropImg6(img_list):
    n = len(img_list)

    assert(n == 6)

    # y_top1 = 630
    # y_bot1 = 1720
    # y_top2 = 2800
    # y_bot2 = 4175

    y_top1 = 550
    y_bot1 = 1850
    y_top2 = 2770
    y_bot2 = 4220

    img1_cropped1 = img_list[0].crop((530, y_top1, 1900, y_bot1))
    img1_cropped2 = img_list[0].crop((530, y_top2, 1900, y_bot2))

    img2_cropped1 = img_list[1].crop((530, y_top1, 1900, y_bot1))
    img2_cropped2 = img_list[1].crop((530, y_top2, 1900, y_bot2))

    img3_cropped1 = img_list[2].crop((530, y_top1, 1900, y_bot1))
    img3_cropped2 = img_list[2].crop((530, y_top2, 1900, y_bot2))

    img4_cropped1 = img_list[3].crop((530, y_top1, 1900, y_bot1))
    img4_cropped2 = img_list[3].crop((530, y_top2, 1900, y_bot2))

    img5_cropped1 = img_list[4].crop((530, y_top1, 1900, y_bot1))
    img5_cropped2 = img_list[4].crop((530, y_top2, 1900, y_bot2))

    img6_cropped1 = img_list[5].crop((530, y_top1, 1900, y_bot1))
    img6_cropped2 = img_list[5].crop((530, y_top2, 1900, y_bot2))


    img1_comb = get_concat_v(img1_cropped1, img1_cropped2)
    img2_comb = get_concat_v(img2_cropped1, img2_cropped2)
    img3_comb = get_concat_v(img3_cropped1, img3_cropped2)
    img4_comb = get_concat_v(img4_cropped1, img4_cropped2)
    img5_comb = get_concat_v(img5_cropped1, img5_cropped2)
    img6_comb = get_concat_v(img6_cropped1, img6_cropped2)

    img12 = get_concat_h(img1_comb, img2_comb)
    img34 = get_concat_h(img3_comb, img4_comb)
    img56 = get_concat_h(img5_comb, img6_comb)

    img1234 = get_concat_h(img12, img34)

    final_img = get_concat_h(img1234, img56)
    return final_img

def cropImg2(img_list):
    n = len(img_list)

    assert(n == 2)

    # y_top1 = 630
    # y_bot1 = 1720
    # y_top2 = 2800
    # y_bot2 = 4175
    
    y_top1 = 550
    y_bot1 = 1850
    y_top2 = 2770
    y_bot2 = 4220

    img1_cropped1 = img_list[0].crop((530, y_top1, 1900, y_bot1))
    img1_cropped2 = img_list[0].crop((530, y_top2, 1900, y_bot2))

    img2_cropped1 = img_list[1].crop((530, y_top1, 1900, y_bot1))
    img2_cropped2 = img_list[1].crop((530, y_top2, 1900, y_bot2))

    img1_comb = get_concat_v(img1_cropped1, img1_cropped2)
    img2_comb = get_concat_v(img2_cropped1, img2_cropped2)

    img12 = get_concat_h(img1_comb, img2_comb)

    return img12

def drawboxplot2(list1, list2, x_label, y_label, title, outputDir, outputName):
    # Define data
    data = [list1, list2]

    # Define figure
    fig, ax = plt.subplots()

    # Draw the boxplots with x_label and y_label
    bplot = ax.boxplot(data, notch=True, labels=x_label)
    ax.set_ylabel(y_label)

    # set title
    ax.set_title(title)

    # save figure
    fig.savefig(os.path.join(outputDir, outputName), dpi=400, format='tif')




