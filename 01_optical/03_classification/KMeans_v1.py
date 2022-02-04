#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
from sklearn import cluster
from sklearn.cluster import MiniBatchKMeans
from osgeo import gdal, gdal_array
import argparse
import sys
import os

#----------------------------------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# KMeans Unsupervised classification
Processes an Unsupervised classification using KMeans and a user defined number of clusters.
To be used to assess which classes can be identified within an image before using a supervised classification method.
**************************************************************************
##Tasks:
- Reads image
- Computes KMeans classification
- Saves output
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#==========================================================
#main
#----------------------------------------------------------

if __name__ == "__main__":
    try:
        print('')
        #----------------------------------------------------------------------------------------------------
        # Retrieval of arguments
        #----------------------------------------------------------------------------------------------------
        parser.add_argument('inputImage',help='Input image')
        parser.add_argument('kMeansClusters',type=int,help='Number of clusters wanted')
        parser.add_argument('outPath',help='out put path of kmeans classification')

        args = parser.parse_args()

        #Pass exceptions and register all drivers
        gdal.UseExceptions()
        gdal.AllRegister()

        #Open input image
        img_ds = gdal.Open(args.inputImage, gdal.GA_ReadOnly)

        #Create empty array with the same dimensions as the input image
        img = np.zeros((img_ds.RasterYSize, img_ds.RasterXSize, img_ds.RasterCount),
                        gdal_array.GDALTypeCodeToNumericTypeCode(img_ds.GetRasterBand(1).DataType))

        for b in range(img.shape[2]):
            img[:, :, b] = img_ds.GetRasterBand(b + 1).ReadAsArray()

        new_shape = (img.shape[0] * img.shape[1], img.shape[2])

        X = img[:, :, :img_ds.RasterCount].reshape(new_shape)

        print('Computing %i clusters'% args.kMeansClusters)
        k_means = MiniBatchKMeans(args.kMeansClusters)
        k_means.fit(X)

        X_cluster = k_means.labels_
        X_cluster = X_cluster.reshape(img[:, :, 0].shape)

        ds = gdal.Open(args.inputImage)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        [cols, rows] = arr.shape

        print('Saving image')
        format = "GTiff"
        driver = gdal.GetDriverByName(format)

        #output image
        outImgBaseName = os.path.basename(args.inputImage)[:-4]
        outKMeans = os.path.join(args.outPath, outImgBaseName +'_KMeans_' + str(args.kMeansClusters) + '_class.tif')

        outDataRaster = driver.Create(outKMeans, rows, cols, 1, gdal.GDT_Byte)
        outDataRaster.SetGeoTransform(ds.GetGeoTransform())
        outDataRaster.SetProjection(ds.GetProjection())

        outDataRaster.GetRasterBand(1).WriteArray(X_cluster)

    except RuntimeError as msg:
        print("\nERROR: ", msg)
