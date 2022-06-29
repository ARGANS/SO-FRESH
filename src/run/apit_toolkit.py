#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing the Automated Polynya Identification Tool run steps. 
"""

# Package loader # 
import os, sys, glob, gc
import functools, itertools, collections
import numpy as np
import numpy.ma as ma
import gdal, rasterio
from skimage.filters import threshold_otsu, threshold_minimum, sobel
from skimage.morphology import binary_erosion, binary_dilation, remove_small_holes, remove_small_objects
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
            if filepath == None: raise RuntimeError(f"A fused product containing all inputted products does not exist. Consider processing this using the 'fusion' step.\nInputted products: {self.products}")
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

        products = dict([(p, []) for p in self.products])
        indices = dict([(p, []) for p in self.products])
        src = rasterio.open(img)
        desc = list(src.descriptions)
        src.close()
        if len(self.products) > 1:  
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

        elif len(self.products) == 1:
            if self.products[0] == "MYD09GA":
                products[self.products[0]].append("None")
                products[self.products[0]].append("None")
                products[self.products[0]].append("None")
                indices[self.products[0]].append(1)
                indices[self.products[0]].append(2)
                indices[self.products[0]].append(3)
            elif self.products[0] == "MYDTBGA":
                products[self.products[0]].append("None")
                indices[self.products[0]].append(1)
            elif self.products[0] == "SIC":
                products[self.products[0]].append("None")
                indices[self.products[0]].append(1)
            elif self.products[0] == "SSS":
                products[self.products[0]].append("None")
                indices[self.products[0]].append(1)
            
        return(products, indices)   

    def expression_builder(self, img, bands_dict, idx_dict):

        """ Build calculation expression ready for gdal_calc. """

        calculation, numbers, letters = [], [], []
        src = rasterio.open(img)
        for bd, ind in zip(bands_dict.keys(), idx_dict.keys()):
            if bd == "MYD09GA" and ind == "MYD09GA":
                # Execute adaptive thresholding.
                array = src.read(idx_dict[ind])
                OT=optical_thresholding(array)
                del array
                opt_threshold=OT.otsu()
                numb = idx_dict[ind]
                let = [chr(ord("@")+x) for x in numb]
                calc = "*".join([f"({b}<={opt_threshold})" for b in let])
                numbers.append(numb)
                letters.append(let)
                calculation.append(calc)
                gc.collect()
            elif bd == "MYDTBGA" and ind == "MYDTBGA":
                numb = idx_dict[ind]
                let = [chr(ord("@")+x) for x in numb]
                calc = "*".join([f"({b}>=265)" for b in let])
                numbers.append(numb)
                letters.append(let)
                calculation.append(calc)
            elif bd == "SIC" and ind == "SIC":
                numb = idx_dict[ind]
                let = [chr(ord("@")+x) for x in numb]
                calc = "*".join([f"({b}<=60)*({b}>=0)" for b in let])
                numbers.append(numb)
                letters.append(let)
                calculation.append(calc)

            elif bd == "SSS" and ind == "SSS":
                print("No SSS data implemented yet.")
        del src
        gc.collect()
        numbers = functools.reduce(lambda x,y:x+y,(numbers))
        letters = functools.reduce(lambda x,y:x+y,(letters))
        band_inputs = " ".join([f"-{l} {img} --{l}_band={n}" for n, l in zip(numbers, letters)])
        calculation = "*".join(calculation)
        return(band_inputs, calculation)

    def mask_ice_sheet_edges(self, o_input, img, index, outfile):

        """ Apply edge detection for ice-sheet boundary detection for mask generation. """
        
        src = rasterio.open(o_input)
        sic = src.read(index)

        mask_src = rasterio.open(img)
        mask_arr = mask_src.read(1)
        # Pull image profile for saving.
        profile = mask_src.profile
        # Extract ice-sheet edge. 1-sobel filter, 2-fill holes in ice-sheet data, 3-run another sobel filter for extracting ice-sheet edge.      
        sobel_=binary_erosion(sobel(sic))
        clean=remove_small_objects(remove_small_holes(sobel_, area_threshold=10000), min_size=500, connectivity=50)
        is_edge=sobel(clean, mode="reflect")
        # Dilate imagery - one dilation is one row of pixels. 
        buffer=binary_dilation(is_edge) #10 km
        buffer=binary_dilation(is_edge) #20 km
        buffer=binary_dilation(is_edge) #30 km
        buffer=binary_dilation(is_edge) #40 km
        buffer=binary_dilation(is_edge) #50 km
        # Mask the output
        foutput = ma.array(mask_arr, mask=buffer, fill_value=0)
        foutput = foutput.filled(fill_value=0)

        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(foutput.astype(rasterio.float32), 1)
        '''      
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 1)
        axs[0].imshow(mask_arr, interpolation='nearest')
        axs[1].imshow(buffer, interpolation='nearest')
        axs[2].imshow(foutput, interpolation='nearest')
        plt.show()
        sys.exit()
        '''

    def apit(self):

        """ Execute the Automated Polynya Identification Tool. """

        process = self.process_organiser()
        outfiles = []
        print("Reading images ...\n")
        for img in tqdm(self.data):
            # Read from inputs the required bands for the tool execution.
            bands_dict, idx_dict = self.image_reader(img)
            # Build an expression to input in to gdal_calc.
            band_inputs, expression = self.expression_builder(img, bands_dict, idx_dict)
            for d in self.dates:
                if d in img:date = d
            outdir = f"{self.data_directory}03_APIT/{'_'.join(self.products)}/{date}/"
            if not os.path.isdir(outdir): os.makedirs(outdir)
            outfile = f"{outdir}03_ESA_SOFRESH_APIT_{'_'.join(self.products)}_{date.replace('/','')}_{self.aoi}.tif"
            # Execute the APIT tool (threshold classification).
            cmd = ("gdal_calc.py --quiet %s --outfile=%s --calc='%s' --overwrite"%(band_inputs, outfile, expression))
            os.system("gdal_calc.py --quiet %s --outfile=%s --calc='%s' --overwrite"%(band_inputs, outfile, expression))
            # If one of the input products is sea-ice concentration, mask out areas on the ice-sheet edge which may have been misidentified. 
            if "SIC" in bands_dict.keys() and idx_dict.keys():
                mask_edge = self.mask_ice_sheet_edges(img, outfile, idx_dict["SIC"][0], outfile)



    