#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk

import os, sys
from datetime import datetime, timedelta
import glob, gdal, functools
import numpy as np
from pathlib import Path
from itertools import groupby
from skimage.filters import threshold_otsu, threshold_minimum
from tqdm import tqdm

class data_selector():
    def __init__(self, data_folder, sdate, edate, products):
        self.data_folder = data_folder
        self.sdate = sdate
        self.edate = edate
        self.products = products


    def collector(self):

        """ Accumulate all data that fits within the date and product boundaries. """

        products = self.product_selector()
        dates = self.dates_builder()
        filepaths = []
        if type(products) == str:
            for d in dates:
                p = glob.glob(products+d+"/02_*.tif")
                if len(p) == 1:filepaths.append(p[0])
                else: print(f"Multiple files are in {products+d}, please leave only the one to be used.")

        ### ADD ELIF STATEMENT AS TO WHETHER THIS PULLS JUST OPTICAL OR THERMAL OR SIC MOSAIC'D PRODUCTS ###

        return(filepaths)


    def dates_builder(self):

        """ Pulls all dates from inputs start and end range. """

        sdate=datetime.strptime(os.path.join(self.sdate), "%Y/%m/%d").date()
        edate=datetime.strptime(os.path.join(self.edate), "%Y/%m/%d").date()
        date_dt=[sdate+timedelta(days=x) for x in range((edate-sdate).days+1)]

        return(["/".join((str(d.year), str("%02d" %d.month), str("%02d" %d.day))) for d in date_dt])


    def product_selector(self):

        """ Selects all products based on entries. """

        try:
            fp = self.data_folder+"fusion/"+"_".join(self.products)+"/"
            if os.path.isdir(fp):return(fp)
        except:
            print("No fused products with that order found")
        else:
            product_fp = []
            for p in self.products:
                if p == "MYD09GA":
                    fp = self.data_folder+"MODIS/"+p+"_006/02_mosaic/"
                    product_fp.append(fp)
                elif p == "MYDTBGA":
                    fp = self.data_folder+"MODIS/"+p+"_006/02_mosaic/"
                    product_fp.append(fp)
                elif p == "SIC":
                    fp = self.data_folder+"AMSR2/sic_extracted/"
                    product_fp.append(fp)
            return(product_fp)



class execute_APIT():
    def __init__(self, data_lst, products):
        self.data_lst=data_lst
        self.products=products

    def band_selector(self):

        """ Based on order entry, find out which bands correspond to which data. """

        products = dict.fromkeys(self.products)
        count=0
        for p in self.products:
            if p == "MYD09GA":
                start=count
                count=count+3
                products[p]=start,count
            elif p == "MYDTBGA" or "SIC":
                start=count
                count=count+1
                products[p]=(start,count)
        return(products)


    def image_reader(self, img, dic):

        """ Reads image and based on product entries, slices data for appropriate band selection """

        ds = gdal.Open(img)
        array = ds.ReadAsArray()
        # Slice array to pull bands.
        sliced_arrays = []
        for key, val in dic.items():
            if key == "MYD09GA":
                optical=array[val[0]:val[1]]
                sliced_arrays.append(optical)
            elif key == "MYDTBGA":
                thermal=array[val[0]:val[1]]
                sliced_arrays.append(thermal)
            elif key == "SIC":
                sic=array[val[0]:val[1]]
                sliced_arrays.append(sic)
        return(sliced_arrays)


    def expression_builder(self, img, reader, prods):

        """ Builds expression to be passed for the threshold calculation. """

        data,numbers,letters, calculation=[],[],[],[]
        for p, d in zip(self.products, reader):
            if p == "MYD09GA":
                stats=optical_thresholding(d)
                opt_threshold=optical_thresholding(d).otsu()
                MYD09GA = d
                MYD09GA_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                MYD09GA_letters = [chr(ord("@")+x) for x in MYD09GA_numbers]
                MYD09GA_calc = "*".join([f"({b}<={opt_threshold})" for b in MYD09GA_letters])
                data.append(MYD09GA)
                numbers.append(MYD09GA_numbers)
                letters.append(MYD09GA_letters)
                calculation.append(MYD09GA_calc)
            elif p == "MYDTBGA":
                MYDTBGA = d
                MYDTBGA_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                MYDTBGA_letters = [chr(ord("@")+x) for x in MYDTBGA_numbers]
                MYDTBGA_calc = "*".join([f"({b}>=265)" for b in MYDTBGA_letters])
                data.append(MYDTBGA)
                numbers.append(MYDTBGA_numbers)
                letters.append(MYDTBGA_letters)
                calculation.append(MYDTBGA_calc)
            elif p == "SIC":
                SIC = d
                SIC_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                SIC_letters = [chr(ord("@")+x) for x in SIC_numbers]
                SIC_calc = "*".join([f"({b}<=60)*({b}>=0)" for b in SIC_letters])
                data.append(SIC)
                numbers.append(SIC_numbers)
                letters.append(SIC_letters)
                calculation.append(SIC_calc)
        numbers = functools.reduce(lambda x,y:x+y,(numbers))
        letters = functools.reduce(lambda x,y:x+y,(letters))
        band_expression = " ".join([f"-{l} {img} --{l}_band={n}" for n, l in zip(numbers, letters)])
        full_expression = "*".join(calculation)

        return(band_expression, full_expression)


    def classification(self, data_folder):

        """ Execute APIT classification """

        outfiles=[]
        prods=self.band_selector()
        for img in tqdm(self.data_lst):
            reader=self.image_reader(img, prods)
            print(reader[0])
            print("----")
            print(reader[1])
            sys.exit()
            if len(self.products) == len(reader):
                band_exp, calc_exp = self.expression_builder(img, reader, prods)
                outdir = str(Path(data_folder).parent)+"/03_APIT/"+str(Path(img).parents[3].name)+"/"+str(Path(img).parents[2].name)+"/"+str(Path(img).parents[1].name)+"/"+str(Path(img).parents[0].name)+"/"
                if not os.path.isdir(outdir): os.makedirs(outdir)
                outfile = outdir+"03"+(os.path.basename(img)[2:])
                outfiles.append(outfile)
                os.system("gdal_calc.py --quiet %s --outfile=%s --calc='%s' --overwrite"%(band_exp, outfile, calc_exp))
                
        return(outfiles)


    def generate_netcdf(self, img_lst, data_folder, sdate, edate):

        """ Generate a netcdf output """

        if sdate[:-6] == edate[:-6]:
            if sdate[5:-3] == edate[5:-3]:
                outdir = str(Path(data_folder).parent)+"/03_APIT/netCDF/"+str(Path(img_lst[0]).parents[3].name)+"/"+str(Path(img_lst[0]).parents[2].name)+"/"+str(Path(img_lst[0]).parents[1].name)+"/"
            else:
                outdir = str(Path(data_folder).parent)+"/03_APIT/netCDF/"+str(Path(img_lst[0]).parents[3].name)+"/"+str(Path(img_lst[0]).parents[2].name)+"/"
        else:
            outdir = str(Path(data_folder).parent)+"/03_APIT/netCDF/"+str(Path(img_lst[0]).parents[3].name)+"/"+sdate[:-6]+"_"+edate[:-6]+"/"
        if not os.path.isdir(outdir): os.makedirs(outdir)

        proj_info=[]
        img_arrays=[]
        for img in img_lst:
            ds = gdal.Open(img)
            cols, rows, proj, geom = ds.RasterXSize, ds.RasterYSize, ds.GetProjection(), ds.GetGeoTransform()
            proj_info.append((ds.RasterXSize, ds.RasterYSize, ds.GetProjection(), ds.GetGeoTransform()))
            img_arrays.append(ds.ReadAsArray())
            ds.close()

        # Check projection info for all images match.
        check = groupby(proj_info)
        if next(check, True) and not next(check, False) == True:
            sum_arr = np.sum(np.stack(img_arrays), axis=0)
            sum_arr[sum_arr==0] = np.nan
            tmptif=outdir+"04_SOFRESH_"+str(Path(img_lst[0]).parents[3].name)+"_POLYNYA_"+"".join(sdate.rsplit("/"))+"_"+"".join(edate.rsplit("/"))+".tif"
            self.array_to_img(sum_arr, proj_info[0], tmptif)

        outnc = outdir+"04_SOFRESH_"+str(Path(img_lst[0]).parents[3].name)+"_POLYNYA_"+"".join(sdate.rsplit("/"))+"_"+"".join(edate.rsplit("/"))+".nc"
        os.system("gdal_translate -q -of netCDF %s %s "%(tmptif, outnc))
        os.system("ncatted -O -a long_name,Band1,o,c,CountFrequency %s"%(outnc))


    @staticmethod
    def array_to_img(array, proj, outfile):
        cols, rows, proj, geom = proj[0],proj[1],proj[2],proj[3]
        outdataset=gdal.GetDriverByName("GTiff").Create(outfile, cols, rows, 1, gdal.GDT_Float32)
        outdataset.SetProjection(proj)
        outdataset.SetGeoTransform(geom)
        outband=outdataset.GetRasterBand(1)
        outband.WriteArray(array)



class optical_thresholding():
    def __init__(self, array, default_thresh=180, max_thresh=190, min_thresh=120):
        self.array=array
        self.default_thresh=default_thresh
        self.max_thresh=max_thresh
        self.min_thresh=min_thresh


    def greyscale(self):

        """ Builds greyscale optical image to be passed to the Otsu threshold. """

        return(np.sum(self.array, axis=0)/self.array.shape[0])


    def otsu(self, bins=256):

        """ Execute an Otsu threshold. """

        try:
            otsu_t=threshold_otsu(image=self.greyscale(), nbins=bins)
            if self.min_thresh < otsu_t < self.max_thresh:
                return(otsu_t)
            else:
                return(self.deault_thresh)
        except:
            return(self.default_thresh)


    def min_threshold(self, bins=256):
        try:
            min_t=threshold_minimum(image=self.greyscale(), nbins=bins)
            if self.min_thresh < min_t < self.max_thresh:
                return(min_t)
            else:
                return(self.default_thresh)
        except:
            return(self.default_thresh)
