#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import sys, os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from PIL import Image
import osr
from tqdm import tqdm
from skimage.filters import threshold_otsu, threshold_minimum
from sklearn.neighbors import KernelDensity
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Undertake a basic threshold classification.
**************************************************************************
##Tasks:
-  Threshold classification where the upper limit value is set.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
class thresholding:
    def __init__(self, data, default_threshold=180, max_thresh=200, min_thresh=160):
        self.data = gdal.Open(data).ReadAsArray()
        self.default_threshold = default_threshold
        self.min_thresh = min_thresh
        self.max_thresh = max_thresh
        self.x = None
        self.y = None
    
    def greyscale(self):
        # Returns greyscale of input image (all array bands are summed and divided by total number)
        return np.sum(self.data, axis=0)/(self.data.shape[0]) 

    def otsu_threshold(self, bins):
        # Thresholding using Otsu's method
        try:
            t = threshold_otsu(image=array, nbins=bins)
            if self.min_thresh < t < self.max_thresh:
                return t
            else:
                return self.default_threshold
        except:
            return self.default_threshold

    def min_threshold(self, bins=256):
        try:
            t = threshold_minimum(image=array, nbins=bins)
            if self.min_thresh < t < self.max_thresh:
                return t
            else:
                return self.default_thresh
        except:
            return self.default_thresh

def threshold_classify(input, output, threshold):
    tmpout = output[:-4] + "thresholdtmp.tif"
    os.system("gdal_calc.py --quiet -A %s --allBands=A --outfile=%s --calc='A<=%s' --overwrite"%(input, tmpout, threshold))
    os.system("gdal_translate -q -b 1 %s %s"%(tmpout, output))
    os.remove(tmpout)

def img_viewer(input):
    img = mpimg.imread(input)
    imgplot = plt.imshow(img)
    plt.show()
#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="Input image - reprojected jpg file.").completer = FilesCompleter(allowednames=(".tif"))
        #parser.add_argument("-t", "--threshold", required=True, type=int, help="The upper threshold value as an interger.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        for img in args.input_img:
            init = thresholding(img)
            #arr = init.greyscale()
            print(init)
            print(init.otsu_threshold(bins=128))

            #print(array)
            #print(otsu)
            #mini = init.min_threshold()
            sys.exit()
            #'''

            plt.hist(array)
            plt.axvline(x=(otsu))
            plt.show()
            sys.exit()
            #'''
            print(test.otsu_threshold(bins=256))
        sys.exit()
        print(f"Undertaking threshold classification for {len(args.input_img)} images")
        for img in tqdm(args.input_img):
            threshold_classify(img, os.path.split(img)[0] + "/03" + os.path.basename(os.path.splitext(img)[0])[2:] + "_threshold_" + str(args.threshold) + ".tif", args.threshold)
        print("Process complete - threshold classifications produced")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
