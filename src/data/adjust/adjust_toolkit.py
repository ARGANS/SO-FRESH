#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing the Automated Polynya Identification Tool pre-processing steps. 
"""

# Package loader # 
import os, sys
from difflib import SequenceMatcher
import gdal
import glob, itertools, osr, functools
import gc # garbagecollector
import numpy as np
import pandas as pd 
import scipy
from skimage import morphology

class myd_adjust:
    def __init__(self, img):
        self.img = img
        #self.tile_extent_csv = csv
        #self.epsg = epsg
    
    def extract_geometry(self, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../", "aux_files/modis_sinusoidal_tiles.csv")), epsg=4326):

        """ Pull extents of tile referenced in the file naming. """

        if self.img.endswith(".jpg"):
            if os.path.basename(self.img).rsplit(".")[3][0] == "h" and os.path.basename(self.img).rsplit(".")[3][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[3][1:3], os.path.basename(self.img).rsplit(".")[3][4:6]
            else:
                raise RuntimeError("Unable to find appropriate matching tiles.")
        elif self.img.endswith(".tif"):
            if os.path.basename(self.img).rsplit(".")[2][0] == "h" and os.path.basename(self.img).rsplit(".")[2][3] == "v":
                h, v = os.path.basename(self.img).rsplit(".")[2][1:3], os.path.basename(self.img).rsplit(".")[2][4:6]
        df = pd.read_csv(csv)
        df_row = df[(df["iv"] == int(v)) & (df["ih"] == int(h))]
        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            ny, nx = img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            ny, nx = img_array.shape[0], img_array.shape[1]
        xmin, ymin, xmax, ymax = float(df_row["lon_min"]), float(df_row["lat_min"]), float(df_row["lon_max"]), float(df_row["lat_max"])
        xres, yres = (xmax-xmin)/nx, (ymax-ymin)/ny

        del df, df_row, img_array, ny, nx, xmax, ymin
        gc.collect()
        return(xmin, xres, 0, ymax, 0, -yres)

    def assign_geometry(self, epsg=4326):

        """ Assign projection to image based on tile name. """

        img_array = gdal.Open(self.img).ReadAsArray()
        if len(img_array.shape) == 3:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            shp, ny, nx = gdal.Open(self.img).RasterCount, img_array.shape[0], img_array.shape[1]
        outfile = os.path.dirname(self.img)+("/02_"+os.path.basename(self.img)[3:-4]+".tif")
        outdataset = gdal.GetDriverByName("GTiff").Create(outfile, ny, nx, shp, gdal.GDT_Float32)
        outdataset.SetGeoTransform(myd_adjust(self.img).extract_geometry())
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        outdataset.SetProjection(srs.ExportToWkt())
        count=0
        for arr in img_array:
            count+=1
            outband=outdataset.GetRasterBand(count)
            outband.WriteArray(arr)
            outband=None
        
        del img_array, outdataset, srs
        gc.collect()

        return(outfile)
    
    @staticmethod
    def resample(img, xres, yres):

        """ Reseample imagery to a given spatial resolution """    
        
        outfile=(os.path.split(img)[0]+"/02"+os.path.basename(img)[3:])
        os.system("gdal_translate -q -tr %s %s -r bilinear %s %s"%(xres, yres, img, outfile))

        return(outfile)

    @staticmethod
    def sort_inputs_by_date(dates, products, output_directories):

        """ Produce dictionary where the keys are each date from input, aiding with efficient download and mosaicing options for each date."""

        mosaic_dir = dict([(d, []) for d in dates])
        mosaic_url = dict([(d, []) for d in dates])

        for prod, dir in zip(products, output_directories):
            # Check matching of #1 product #2 date & #3 tile in name
            p_split, d_split = prod.rsplit("/"), dir.rsplit("/")
            p_match, d_match = SequenceMatcher(None, p_split[4], d_split[5]).ratio(), SequenceMatcher(None, p_split[5], os.path.join(d_split[8],d_split[9],d_split[10])).ratio()
            if d_split[7] in p_split[6] and p_match >= 0.89 and d_match >= 0.79:     
                mosaic_dir[os.path.join(d_split[8],d_split[9],d_split[10])].append(dir)
                mosaic_url[os.path.join(d_split[8],d_split[9],d_split[10])].append(prod)

            else:
                raise RuntimeError("The order of the list is not consistent - files appear to be missing.")

        return(mosaic_dir, mosaic_url)

    @staticmethod
    def mosaic(data_directory, imgs, product, date, aoi):

        if product[0] == "MYD09GA":product = "MYD09GA_061"
        elif product[0] == "MYDTBGA":product = "MYDTBGA_006"
        outdir=(data_directory+"02_data/MODIS/"+product+"/02_mosaic/")+date
        if not os.path.isdir(outdir): os.makedirs(outdir)
        txt=outdir+"/"+product+"_"+("".join(date.rsplit("/")))+"_files4merge.txt"
        with open(txt, "w") as f:
            for item in imgs:
                f.write("%s\n" % item)
        outvrt=outdir+("/02a_"+product+"_"+("".join(date.rsplit("/")))+"_"+aoi+".vrt")
        os.system("gdalbuildvrt -q -overwrite --config GDAL_VRT_TRUSTED_MODULES YES %s -srcnodata 0 -input_file_list %s"%(outvrt, txt))
        # Modify the xml (VRT) file for claculating mean pixel value where images overlap.
        vrtimg = gdal.Open(outvrt)
        vrtcount = vrtimg.RasterCount
        myd_adjust.vrt_mean_overlap(outvrt, vrtcount)
        outfile=outdir+("/02_"+product+"_"+("".join(date.rsplit("/")))+"_"+aoi+".tif")
        os.system("gdal_translate -q --config GDAL_VRT_ENABLE_PYTHON YES %s %s"%(outvrt, outfile))

        os.remove(txt)
        del vrtimg
        gc.collect()
        return(outfile)


    @staticmethod
    def vrt_mean_overlap(vrt, band):
        if band == 1: # If image is thermal.
            header = """  <VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">"""
            contents = """
            <PixelFunctionType>{0}</PixelFunctionType>
            <PixelFunctionLanguage>Python</PixelFunctionLanguage>
            <PixelFunctionCode><![CDATA[{1}]]>
            </PixelFunctionCode>"""
            lines = open(vrt, "r").readlines()
            lines[3] = header
            lines.insert(4, contents.format("average", """
import numpy as np
import sys
def average(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize,raster_ysize, buf_radius, gt, **kwargs):
    div = np.zeros(in_ar[0].shape)
    
    for i in range(len(in_ar)):
        div += (in_ar[i] != 0)
    div[div == 0] = 1
    
    y = np.sum(in_ar, axis = 0, dtype = 'Float32')
    y = y / div
    
    np.clip(y,0,300, out = out_ar)"""))
            open(vrt, 'w').write("".join(lines))

        elif band == 3: # If image is RGB.
            for b in range(1,4):
                header_to_find = f"""  <VRTRasterBand dataType="Float32" band="{b}">"""
                header = f"""  <VRTRasterBand dataType="Float32" band="{b}" subClass="VRTDerivedRasterBand">"""
                contents = """
                <PixelFunctionType>{0}</PixelFunctionType>
                <PixelFunctionLanguage>Python</PixelFunctionLanguage>
                <PixelFunctionCode><![CDATA[{1}]]>
                </PixelFunctionCode>"""
                lines = open(vrt, "r").readlines()
                count=0
                for l in lines:
                    if header_to_find in l:
                        lines[count] = header
                        lines.insert(count+1, contents.format("average", """
import numpy as np
import sys
def average(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize,raster_ysize, buf_radius, gt, **kwargs):
    div = np.zeros(in_ar[0].shape)
    
    for i in range(len(in_ar)):
        div += (in_ar[i] != 0)
    div[div == 0] = 1
    
    y = np.sum(in_ar, axis = 0, dtype = 'Float32') #total of all overlapping pixels
    y = y / div # total of pixel, divided by number of layers
    
    np.clip(y,0,255, out = out_ar)"""))
                        open(vrt, 'w').write("".join(lines))
                    else:
                        pass
                    count = count + 1

    @staticmethod
    def aoi_tile_identifier(hmin, hmax, vmin, vmax, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../", "aux_files/modis_sinusoidal_tiles.csv"))):
                
        """ Identify all tile posibilities from entries. """
        
        h = [hmin+i for i in range((hmax - hmin)+1)]
        v = [vmin+i for i in range((vmax - vmin)+1)]
        df = pd.read_csv(csv)
        posibilities = list(itertools.product(h, v))
        tiles=[]
        for p in posibilities:
            df_row = df[(df["iv"] == int(p[1])) & (df["ih"] == int(p[0]))]
            if (df_row["lon_min"].values != int(-999)) and (df_row["lon_max"].values != int(-999)) and (df_row["lat_min"].values != int(-99)) and (df_row["lat_max"].values != int(-99)):
                tile = "".join(("h",str(p[0]),"v",str(p[1])))
                if tile == "h16v15": # This tile is missing from the MODIS records. 
                    continue
                tiles.append(tile)
        return(tiles)

class myd09ga_adjust:
    def __init__(self, img):
        self.img = img

    def reproject(self):
        
        outfile = myd_adjust(self.img).assign_geometry()

        return(outfile)

class mydtbga_adjust:
    def __init__(self, img):
        self.img = img

class amsr2_adjust:
    def __init__(self, img, resample, epsg=4326):
        self.img = img
        self.resample = resample
        self.epsg = epsg

    def extract_from_netcdf(self):

        """ Pull the 'ice_conc' variable from the netCDF file and convert it to a tif. """

        date = self.img.rsplit("-")[-2]
        year, month, day = date[:-4], date[4:-2], date[6:]
        outdir = os.path.join(os.path.abspath(os.path.join(os.path.dirname(self.img), "../../../02_sic")),year,month,day)
        if not os.path.isdir(outdir): os.makedirs(outdir)
        os.system("gdal_translate -q -ot Float32 NETCDF:%s:ice_conc %s"%(self.img, outdir+"/01a_"+os.path.basename(os.path.splitext(self.img)[0])[3:]+".tif"))

        return(outdir+"/01a_"+os.path.basename(os.path.splitext(self.img)[0])[3:]+".tif")

    def reproject(self, img):

        """ Assign known CRS projection to Sea Ice Concentration file. """

        img_open = gdal.Open(img)
        xsize=img_open.RasterXSize
        ysize=img_open.RasterYSize
        meta = img_open.GetMetadata()
        del img_open
        lat_max=meta["NC_GLOBAL#geospatial_lat_max"]
        lat_min=meta["NC_GLOBAL#geospatial_lat_min"]
        lon_max=meta["NC_GLOBAL#geospatial_lon_max"]
        lon_min=meta["NC_GLOBAL#geospatial_lon_min"]
        #if bool(self.resample) == True:outfile=os.path.dirname(img)+"/02a_"+os.path.basename(os.path.splitext(img)[0])[4:]+".tif"
        outfile=os.path.dirname(img)+"/02a_"+os.path.basename(os.path.splitext(img)[0])[4:]+".tif"
        os.system("gdalwarp -q -overwrite -ot Float32 -t_srs EPSG:%s -te %s %s %s %s -ts %s %s %s %s"%(self.epsg, lon_min, lat_min, lon_max, lat_max, xsize, ysize, img, outfile))
        
        # Clear system.
        #os.remove(img) # Remove the 01a file.
        gc.collect()
        
        return(outfile)

    @staticmethod
    def resample(img, xres, yres):

        """ Resample AMSR2 product."""
        
        sf = float(gdal.Open(img).GetMetadata()["ice_conc#scale_factor"])
        outvrt=os.path.dirname(img)+"/02b_"+os.path.basename(os.path.splitext(img)[0])[4:]+".vrt"
        outfile=os.path.dirname(img)+"/02c_"+os.path.basename(os.path.splitext(img)[0])[4:]+".tif"
        #xres, yres = self.resample[0], self.resample[1]
        os.system("gdalbuildvrt -q -tr %s %s -r bilinear %s %s"%(xres, yres, outvrt, img))
        os.system("gdal_calc.py --quiet -A %s --type='Float32' --outfile %s --calc 'A*%s'"%(outvrt, outfile, sf))
        
        # Clear system.
        os.remove(img)
        os.remove(outvrt)
        del sf
        gc.collect()
        return(outfile)

    @staticmethod
    def mask(img):

        """ Generate mask to remove all data which is not the ice-sheet. """

        img_open = gdal.Open(img)
        cols, rows, proj, geom = img_open.RasterXSize, img_open.RasterYSize, img_open.GetProjection(), img_open.GetGeoTransform()
        arr = img_open.ReadAsArray()
        if os.path.basename(img).startswith("02c"):
            # Resampled image.
            arr[np.where(arr>100)]=0 # Anything beyond 100, set to 0
        elif os.path.basename(img).startswith("02a"):
            # Not resampled image.
            arr = arr/100 # Change magnitude so values are equal to 0
            arr[np.where(arr>100)]=0 # Anything beyond 100, set to 0
            arr[np.where(arr<0)]=-9999 #-32767.0 # Where this value was divided by 100, set back to the 'no data value'.

        mask_arr = morphology.remove_small_objects(scipy.ndimage.binary_fill_holes(arr), min_size=500)
        arr[np.where(mask_arr<=False)]=-9999
        outfile = os.path.dirname(img)+"/02_"+os.path.basename(os.path.splitext(img)[0])[4:]+".tif"
        outdataset=gdal.GetDriverByName("GTiff").Create(outfile, cols, rows, 1, gdal.GDT_Float32)
        outdataset.SetProjection(proj)
        outdataset.SetGeoTransform(geom)
        outband=outdataset.GetRasterBand(1)
        outband.WriteArray(arr)
        del img_open, arr
        os.remove(img)
        gc.collect()