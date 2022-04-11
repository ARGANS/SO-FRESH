#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing functions for pre-processing MODIS and AMSR2 products.
"""
# packages
import os, sys
import gdal, glob, itertools, osr, functools
import numpy as np
import pandas as pd
import scipy
import skimage
from skimage import morphology
from datetime import datetime, timedelta

class aux_func():
    @staticmethod
    def rename(path):

        """ Rename MODIS imagery by adding "01" at the begining for easy distinguishability. """

        for root, dirs, files in os.walk(path):
            for i in files:
                if i.startswith("BROWSE") and i.endswith(".jpg"):
                    os.rename(os.path.join(root, i), os.path.join(root, "01_"+i))
                elif i.startswith("MYDTBGA") and i.endswith(".hdf"):
                    os.rename(os.path.join(root, i), os.path.join(root, "01_"+i))
                elif i.startswith("ESACCI") and i.endswith(".nc"):
                    os.rename(os.path.join(root, i), os.path.join(root, "01_"+i))

    @staticmethod
    def array_to_tiff(array, image, outfile):

        """ Write a 2D array to a 1 band image """

        img=gdal.Open(image)
        cols, rows, proj, geom = img.RasterXSize, img.RasterYSize, img.GetProjection(), img.GetGeoTransform()
        outdataset=gdal.GetDriverByName("GTiff").Create(outfile, cols, rows, 1, gdal.GDT_Float32)
        outdataset.SetProjection(proj)
        outdataset.SetGeoTransform(geom)
        outband=outdataset.GetRasterBand(1)
        outband.WriteArray(array)

class MYD09GA_preprocess():
    def __init__(self, img, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "aux_files/modis_sinusoidal_tiles.csv")), epsg=4326):
        self.img = img
        self.tile_extent_csv = csv
        self.epsg = epsg

    def extract_geometry(self):

        """ Pull extents of tile referenced in the file naming. """

        if self.img.endswith(".jpg"):
            if os.path.basename(self.img).rsplit(".")[3][0] == "h" and os.path.basename(self.img).rsplit(".")[3][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[3][1:3], os.path.basename(self.img).rsplit(".")[3][4:6]
            else:
                raise RuntimeError("Unable to find appropriate matching tiles.")
        elif self.img.endswith(".tif"):
            if os.path.basename(self.img).rsplit(".")[2][0] == "h" and os.path.basename(self.img).rsplit(".")[2][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[2][1:3], os.path.basename(self.img).rsplit(".")[2][4:6]
        df = pd.read_csv(self.tile_extent_csv)
        df_row = df[(df["iv"] == int(v)) & (df["ih"] == int(h))]
        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            ny, nx = img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            ny, nx = img_array.shape[0], img_array.shape[1]
        xmin, ymin, xmax, ymax = float(df_row["lon_min"]), float(df_row["lat_min"]), float(df_row["lon_max"]), float(df_row["lat_max"])
        xres, yres = (xmax-xmin)/nx, (ymax-ymin)/ny

        return(xmin, xres, 0, ymax, 0, -yres)

    def assign_geometry(self):

        """ Assign projection to image based on tile name. """

        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[0], img_array.shape[1]
        outfile = os.path.dirname(self.img)+("/02_"+os.path.basename(self.img)[3:-4]+".tif")
        outdataset = gdal.GetDriverByName("GTiff").Create(outfile, ny, nx, shp, gdal.GDT_Float32)
        outdataset.SetGeoTransform(modis_preprocess(self.img).extract_geometry())
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.epsg)
        outdataset.SetProjection(srs.ExportToWkt())
        count=0
        for arr in img_array:
            count+=1
            outband=outdataset.GetRasterBand(count)
            outband.WriteArray(arr)
            outband=None

class MYDTBGA_preprocess():
    def __init__(self, img, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "aux_files/modis_sinusoidal_tiles.csv")), epsg=4326):
        self.img = img
        self.tile_extent_csv = csv
        self.epsg = epsg

    def extract_geometry(self):

        """ Pull extents of tile referenced in the file naming. """

        if self.img.endswith(".jpg"):
            if os.path.basename(self.img).rsplit(".")[3][0] == "h" and os.path.basename(self.img).rsplit(".")[3][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[3][1:3], os.path.basename(self.img).rsplit(".")[3][4:6]
            else:
                raise RuntimeError("Unable to find appropriate matching tiles.")
        elif self.img.endswith(".tif"):
            if os.path.basename(self.img).rsplit(".")[2][0] == "h" and os.path.basename(self.img).rsplit(".")[2][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[2][1:3], os.path.basename(self.img).rsplit(".")[2][4:6]
        df = pd.read_csv(self.tile_extent_csv)
        df_row = df[(df["iv"] == int(v)) & (df["ih"] == int(h))]
        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            ny, nx = img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            ny, nx = img_array.shape[0], img_array.shape[1]
        xmin, ymin, xmax, ymax = float(df_row["lon_min"]), float(df_row["lat_min"]), float(df_row["lon_max"]), float(df_row["lat_max"])
        xres, yres = (xmax-xmin)/nx, (ymax-ymin)/ny

        return(xmin, xres, 0, ymax, 0, -yres)

    def assign_geometry(self):

        """ Assign projection to image based on tile name. """

        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[0], img_array.shape[1]
        outfile = os.path.dirname(self.img)+("/02_"+os.path.basename(self.img)[3:-4]+".tif")
        outdataset = gdal.GetDriverByName("GTiff").Create(outfile, ny, nx, shp, gdal.GDT_Float32)
        outdataset.SetGeoTransform(modis_preprocess(self.img).extract_geometry())
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.epsg)
        outdataset.SetProjection(srs.ExportToWkt())
        count=0
        for arr in img_array:
            count+=1
            outband=outdataset.GetRasterBand(count)
            outband.WriteArray(arr)
            outband=None

    def extract_MYDTBGA(self):

        """ Pull thermal band 4 from MYDTBGA HDF product and convert to tif (MODIS Band 32). """

        process_dir = os.path.split(os.path.abspath(self.img))[0]+("/01a_"+os.path.basename(os.path.splitext(self.img)[0])[3:])
        if not os.path.isdir(process_dir): os.mkdir(process_dir)
        os.system("gdal_translate -q -f GTIFF -sds %s %s"%(self.img, (process_dir+"/01b_"+os.path.basename(os.path.splitext(self.img)[0])[3:]+".tif")))

        return(process_dir+"/01b_"+os.path.basename(os.path.splitext(self.img)[0])[3:])

    def normalise(self):

        """ Normalise thermal image and assign projection. """

        outfile = os.path.dirname(os.path.dirname(self.img))+("/02_"+os.path.basename(self.img)[4:-6]+".tif")
        o_img = gdal.Open(self.img)
        arr_img = o_img.ReadAsArray()
        shp, nx, ny = o_img.RasterCount, arr_img.shape[0], arr_img.shape[1]
        outdataset = gdal.GetDriverByName("GTiff").Create(outfile, ny, nx, shp, gdal.GDT_Float32)
        outdataset.SetGeoTransform(modis_preprocess(self.img).extract_geometry())
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.epsg)
        outdataset.SetProjection(srs.ExportToWkt())
        outband = outdataset.GetRasterBand(1)
        outband.SetDescription("Thermal")
        outband.WriteArray(arr_img*0.01)


class MODIS():
    def __init__(self, data_fp, product):
        self.data_fp=data_fp
        self.product=product

    def build_filepath(self, sdate, edate):

        """ Find files based on dates & products for mosaic building. """

        sdate=datetime.strptime(os.path.join(sdate), "%Y/%m/%d").date()
        edate=datetime.strptime(os.path.join(edate), "%Y/%m/%d").date()
        dates=[sdate+timedelta(days=x) for x in range((edate-sdate).days+1)]
        if self.product == "MYD09GA" or "MYDTBGA":
            ####### VERSION TO BE REMOVED ONCE 061 EXISTS!!! ###########
            ### Create function where used can specify "antarctica" and it will select tiles for that region
            ### And they can specify specific tiles to merge
            tiles = sorted(glob.glob((self.data_fp+"MODIS/"+self.product+"_006/"+"01_tiles"+"/*")))
            dates = ["/".join((str(d.year), str("%02d" %d.month), str("%02d" %d.day))) for d in dates]
        else: raiseRuntimeError("A MODIS product was not picked up, please check product entry")
        mosaic=[]
        for d in dates:
            imgs=[]
            for t in tiles:
                files_lst = glob.glob((t+"/"+d+"/02_*.tif"))
                if len(files_lst) == 1: path=files_lst[0]
                imgs.append(path)
            mosaic.append(imgs)
        return(mosaic, dates)

    def build_mosaic(self, imgs, date, r):
        if self.product == "MYD09GA" or "MYDTBGA":
            #### VERSION TO BE REMOVED ONCE 061 exists for MYDTBGA #####
            outdir=(self.data_fp+"MODIS/"+self.product+"_006/02_mosaic/")+date
            if not os.path.isdir(outdir): os.makedirs(outdir)
            if r == True:outfile=outdir+("/02a_"+self.product+"_"+("".join(date.rsplit("/")))+"_ANTARCTICA.tif")
            else:outfile=outdir+("/02_"+self.product+"_"+("".join(date.rsplit("/")))+"_ANTARCTICA.tif")
            os.system("gdal_merge.py -q -o %s %s"%(outfile, " ".join(imgs)))

            return(outfile)

    def resample(self, img):
        xres, yres = 0.1, 0.1
        outfile=(os.path.split(img)[0]+"/02"+os.path.basename(img)[3:])
        os.system("gdal_translate -q -tr %s %s -r bilinear %s %s"%(xres, yres, img, outfile))


class amsr2_preprocess():
    def __init__(self, img, epsg=4326):
        self.img = img
        self.epsg = epsg

    def extract_from_netcdf(self):

        """ Pull the 'ice_conc' variable from the netCDF file and convert it to a tif. """

        date = self.img.rsplit("-")[-2]
        year, month, day = date[:-4], date[4:-2], date[6:]
        outdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname(self.img), "../../../sic_extracted")),year,month,day)
        if not os.path.isdir(outdir): os.mkdir(outdir)
        os.system("gdal_translate -q -ot Float32 NETCDF:%s:ice_conc %s"%(self.img, outdir+"/01a_"+os.path.basename(os.path.splitext(self.img)[0])[3:]+".tif"))

        return(outdir+"/01a_"+os.path.basename(os.path.splitext(self.img)[0])[3:]+".tif")

    def reproject(self, r):

        """ Assign known CRS projection to Sea Ice Concentration file. """

        img = gdal.Open(self.img)
        meta = img.GetMetadata()
        lat_max=meta["NC_GLOBAL#geospatial_lat_max"]
        lat_min=meta["NC_GLOBAL#geospatial_lat_min"]
        lon_max=meta["NC_GLOBAL#geospatial_lon_max"]
        lon_min=meta["NC_GLOBAL#geospatial_lon_min"]
        xsize=img.RasterXSize
        ysize=img.RasterYSize
        if r == True:outfile=os.path.dirname(self.img)+"/02a_"+os.path.basename(os.path.splitext(self.img)[0])[4:]+".tif"
        else:outfile=os.path.dirname(self.img)+"/02_"+os.path.basename(os.path.splitext(self.img)[0])[4:]+".tif"
        os.system("gdalwarp -q -overwrite -ot Float32 -t_srs EPSG:%s -te %s %s %s %s -ts %s %s %s %s"%(self.epsg, lon_min, lat_min, lon_max, lat_max, xsize, ysize, self.img, outfile))

        return(outfile)

    def resample(self):

        """ Resample AMSR2 product."""

        xres, yres = 0.1, 0.1
        sf = float(gdal.Open(self.img).GetMetadata()["ice_conc#scale_factor"])
        outvrt=os.path.dirname(self.img)+"/02b_"+os.path.basename(os.path.splitext(self.img)[0])[4:]+".vrt"
        outfile=os.path.dirname(self.img)+"/02c_"+os.path.basename(os.path.splitext(self.img)[0])[4:]+".tif"
        os.system("gdalbuildvrt -q -tr %s %s -r bilinear %s %s"%(xres, yres, outvrt, self.img))
        os.system("gdal_calc.py --quiet -A %s --type='Float32' --outfile %s --calc 'A*%s'"%(outvrt, outfile, sf))

        return(outfile)

    def mask(self):

        """ Generate mask to remove all data which is not the ice-sheet. """

        arr = gdal.Open(self.img).ReadAsArray()
        arr[np.where(arr>100)]=0
        mask_arr=morphology.remove_small_objects(scipy.ndimage.binary_fill_holes(arr), min_size=500)
        arr[np.where(mask_arr<=False)]=-9999
        aux_func().array_to_tiff(arr, self.img, os.path.dirname(self.img)+"/02_"+os.path.basename(os.path.splitext(self.img)[0])[4:]+".tif")