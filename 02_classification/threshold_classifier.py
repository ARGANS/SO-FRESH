#!/usr/bin/env python

import argparse
import gdal
import sys, os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import osr
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Undertake a basic threshold classification.
**************************************************************************
##Tasks:
-  Threshold classification where the upper limit value is selected.
-
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def threshold_classify(input, output, threshold):
    tmpout = output[:-4] + "thresholdtmp.tif"
    os.system("gdal_calc.py -A %s --allBands=A --outfile=%s --calc='A<=%s' --overwrite"%(input, tmpout, threshold))
    os.system("gdal_translate -b 1 %s %s"%(tmpout, output))
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
        parser.add_argument("-i", "--input-img", required=True, help="Input image to be classified.")
        parser.add_argument("-o", "--output-img", required=True, help="Output image filepath and name.")
        parser.add_argument("-t", "--threshold", required=True, type=int, help="The upper threshold value as an interger.")
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        threshold_classify(args.input_img, args.output_img, args.threshold)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
