#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing the Automated Polynya Identification Tool run steps. 
"""

# Package loader # 
import os, sys, glob, gc
import functools, itertools, collections
import numpy as np
import rasterio
from skimage.filters import threshold_otsu, threshold_minimum
from tqdm import tqdm


class apit_preparation():
    def __init__(self, data_directory, dates, products, aoi):
        self.data_directory = data_directory
        self.dates = dates
        self.products = products
        self.aoi = aoi

    def data_collector(self):

        """ Pull all data based on inputs for the execution of APIT. """

        products_4_apit = []
        product_path = self.filepath_builder()
        for d in self.dates:
            products_4_apit.append(glob.glob(f"{product_path+d}/02_*{self.aoi}.tif") + glob.glob(f"{product_path+d}/02_*.tif"))
        return(functools.reduce(lambda x,y:x+y,(products_4_apit)))


    def filepath_builder(self):
        
        """ Identify appropriate filepath for acquiring data for implementation in to APIT. """

        if len(self.products) == 1:
            if self.products[0] == "MYD09GA":
                filepath = self.data_directory+f"02_data/MODIS/{self.products[0]}_061/02_mosaic/"
            elif self.products[0] == "MYDTBGA":
                filepath = self.data_directory+f"02_data/MODIS/{self.products[0]}_006/02_mosaic/"
            elif self.products[0] == "SIC":
                filepath = self.data_directory+f"02_data/AMSR2/02_{self.products[0].lower()}/"
        elif len(self.products) > 1:
            # Search for fusion product which contains all requested products.
            filepath = None
            for fp in glob.glob(self.data_directory+f"02_data/fusion/*"):
                if all(p in os.path.split(fp)[1] for p in self.products):
                    filepath = fp+"/"
            if filepath == None: raise RuntimeError("A fused product containing all inputted products does not exist. Consider processing this using the 'fusion' step.")
            if not os.path.isdir(filepath):
                raise RuntimeError(f"The fused data you wish to implement does not exist, either put the products in the correct order or create them using 'fusion'. \n{filepath}")

        return(filepath)

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

class apit_execute():
    def __init__(self, data_directory, dates, products, aoi, data):
        self.data_directory = data_directory
        self.dates = dates
        self.products = products
        self.aoi = aoi
        self.data = data

    def process_organiser(self):

        """ Based on required data, adapt what additional steps may be required based on input data. """

        if collections.Counter(self.products) == collections.Counter(["MYD09GA", "MYDTBGA"]):
            process = "apit and filter"
        elif collections.Counter(self.products) == collections.Counter(["MYD09GA"]):
            process = "apit and filter"
        elif collections.Counter(self.products) == collections.Counter(["MYDTBGA"]):
            process = "apit and filter"
        else:
            process = "apit"

        return(process)

    def image_reader(self, img):

        """ Read input image and pull required bands band on input products. """

        if len(self.products) > 1:  
            products = dict([(p, []) for p in self.products])
            indices = dict([(p, []) for p in self.products])
            arrays = dict([(p, []) for p in self.products])
            src = rasterio.open(img)
            profile = src.profile
            desc = list(src.descriptions)
            src.close()
            # Pull product band names.
            for p in self.products:
                if p == "MYD09GA":
                    products[p].append([s for s in desc if "Red" in s][0])
                    products[p].append([s for s in desc if "Green" in s][0])
                    products[p].append([s for s in desc if "Blue" in s][0])
                elif p == "MYDTBGA":
                    products[p].append([s for s in desc if "Thermal" in s][0])
                elif p == "SIC":
                    products[p].append([s for s in desc if "SIC" in s][0])
                elif p == "SSS":
                    products[p].append([s for s in desc if "SSS" in s][0])
            
            # Pull the index of each band for array extraction.
            for p in products.keys():
                for v in products[p]:
                    indices[p].append(list(i+1 for i, s in enumerate(desc) if v in s)[0])
            '''
            band_out, array_out = [],[]
            for pk, ik in zip(products.keys(), indices.keys()):
                for p, i in zip(products[pk], indices[ik]):
                    array_out.append(src.read(i+1))
                    band_out.append(p)
            '''
        elif len(self.products) == 1:
            print("Setup code to run APIT with one product.")


        return(products, indices, profile)   

    def expression_builder(self, img, bands_dict, idx_dict):


        #https://rasterio.readthedocs.io/en/latest/topics/calc.html

        calculation = []
        src = rasterio.open(img)
        for bd, ind in zip(bands_dict.keys(), idx_dict.keys()):
            if bd == "MYD09GA" and ind == "MYD09GA":
                array = src.read(idx_dict[ind])
                OT=optical_thresholding(array)
                opt_threshold=OT.otsu()
                calc = " ".join([f"(* (<= (read 1 {b}) {opt_threshold}) 1)" for b in idx_dict[ind]])
                calculation.append(calc)
                del array
                gc.collect()
            elif bd == "MYDTBGA" and ind == "MYDTBGA":
                calc = " ".join([f"(* (>= (read 1 {b}) 265) 1)" for b in idx_dict[ind]])
                calculation.append(calc)
            elif bd == "SIC" and ind == "SIC":
                ###### WHY ARE YOU NOT WORKING??? #####
                ## SYNTAX FOR SIC BROKEN??? ###
                #array = src.read(idx_dict[ind])
                #array_masked = np.ma.masked_array(array, mask=(array == -9999))
                #calc = " ".join([f"(* (>= (read 1 {b}) 0) 1) (* (<= (read 1 {b}) 60) 1) (where (== (read 1 {b}) -9999) 0)" for b in idx_dict[ind]])
                #calc = " ".join([f"(* (<= (read 1 {b}) 60) 1)" for b in idx_dict[ind]])
                calc = " ".join([f"(* (<= (read 1 {b}) 60) 1) (>= (read 1 {b}) 0) 1)" for b in idx_dict[ind]])
                #calc = " ".join([f"(where (>= (read 1 {b}) 0) 1) (where (<= (read 1 {b}) 60) 1) (where (== (read 1 {b}) -9999) 0)" for b in idx_dict[ind]])
                #masked = np.ma.masked_where(array == -9999, array)
                #print(type(masked))
                #print(calc)
                #sys.exit()
                calculation.append(calc)
            elif bd == "SSS" and ind == "SSS":
                print("No SSS data implemented yet.")
        del src
        gc.collect()
        
        return(" ".join(calculation))
        #sys.exit()

    def apit(self):
        process = self.process_organiser()

        print("Reading images ...\n")
        for img in tqdm(self.data):
            bands_dict, idx_dict, profile = self.image_reader(img)
            expression = self.expression_builder(img, bands_dict, idx_dict)
            for d in self.dates:
                if d in img:date = d
            outdir = f"{self.data_directory}/03_APIT/{'_'.join(self.products)}/{date}/"
            if not os.path.isdir(outdir): os.makedirs(outdir)
            outfile = f"{outdir}03_ESA_SOFRESH_APIT_{'_'.join(self.products)}_{date.replace('/','')}_{self.aoi}.tif"
            #print(expression)
            #print(img)
            #print(outfile)
            #sys.exit()
            cmd = ('rio calc --overwrite "%s" %s %s'%(expression, img, outfile))
            print(cmd)
            os.system(cmd)
            
            sys.exit()


    