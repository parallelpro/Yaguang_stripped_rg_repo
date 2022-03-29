from matplotlib import rcParams
rcParams["figure.dpi"] = 100
rcParams["savefig.dpi"] = 100

import numpy as np
import corner
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from lightkurve import search_lightcurvefile
import lightkurve as lk
import os 
from matplotlib.colors import ListedColormap
import sys

# define colors
red = '#ca0020'
lightred = '#f4a582'
yellow = '#fee090'
lightblue = '#92c5de'
blue = '#0571b0'
darkgrey = 'darkgray'
black = 'k'


import os 
work_dir = '../'
sys.path.append(work_dir)
overleaf_path = work_dir+'paper/'

fontsize = 7 # minimum fontsize is 5
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.size"] = fontsize #7.5
matplotlib.rcParams["legend.fontsize"] = 2.5#7.5
matplotlib.rcParams['text.usetex'] = True#False #True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}\usepackage{sfmath}'
matplotlib.rcParams['axes.labelsize'] = fontsize#7
matplotlib.rcParams['xtick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.direction']='out'
matplotlib.rcParams['ytick.major.size']=3.0
matplotlib.rcParams['ytick.minor.size']=2.0
matplotlib.rcParams['xtick.direction']='out'
matplotlib.rcParams['xtick.major.size']=3.0
matplotlib.rcParams['xtick.minor.size']=2.0
matplotlib.rcParams['font.family'] = 'sans-serif' #'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica'] #'Helvetica' 


# nature journal size in pt
# https://www.nature.com/documents/NRJs-guide-to-preparing-final-artwork.pdf
columnwidth = 249.449 #88mm
textwidth = 510.236 #180mm

def nature_size(column="one", square=False, ratio=None):
    # Thanks Dan!
    # Parameters:
    # column: "one" or "double"
    # square: True or False
    # ratio: height/width

    inches_per_pt = 1.0/72.00              # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0     # Most aesthetic ratio
    if (ratio == None): ratio = golden_mean
    if (column == "one"):
        fig_width_pt = columnwidth
    elif (column == "double"):
        fig_width_pt = textwidth
    else:
        raise ValueError("column should be one of ``one'' or ``double''. ")
    fig_width = fig_width_pt*inches_per_pt # Figure width in inches
    if square:
        fig_height = fig_width
    else:
        fig_height = fig_width*ratio
    return [fig_width,fig_height]

