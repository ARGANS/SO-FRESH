
  
#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing the Automated Polynya Identification Tool pre-processing steps. 
"""

# Package loader # 
import os, sys

class myd_adjust:
    def __init__():
        print()
    
    def extract_geometry(img, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../", "aux_files/modis_sinusoidal_tiles.csv")), epsg=4326):

        """ Pull extents of tile referenced in the file naming. """

        if img.endswith(".jpg"):
            if os.path.basename(img).rsplit(".")[3][0] == "h" and os.path.basename(img).rsplit(".")[3][3] == "v":
                h, v = os.path.basename(img).rsplit(".")[3][1:3], os.path.basename(img).rsplit(".")[3][4:6]
            else:
                raise RuntimeError("Unable to find appropriate matching tiles.")
        elif img.endswith(".tif"):
            if os.path.basename(img).rsplit(".")[2][0] == "h" and os.path.basename(img).rsplit(".")[2][3] == "v":
                h, v = os.path.basename(img).rsplit(".")[2][1:3], os.path.basename(img).rsplit(".")[2][4:6]
        df = pd.read_csv(csv)
        df_row = df[(df["iv"] == int(v)) & (df["ih"] == int(h))]
        img_array = gdal.Open(img).ReadAsArray()
        if len(img_array.shape) == 3:
            ny, nx = img_array.shape[1], img_array.shape[2]
        elif len(img_array.shape) == 2:
            ny, nx = img_array.shape[0], img_array.shape[1]
        xmin, ymin, xmax, ymax = float(df_row["lon_min"]), float(df_row["lat_min"]), float(df_row["lon_max"]), float(df_row["lat_max"])
        xres, yres = (xmax-xmin)/nx, (ymax-ymin)/ny

        return(xmin, xres, 0, ymax, 0, -yres)

class myd09ga_adjust:
    def __init__(self, img):
        self.img = img
        self.tile_extent_csv = csv
        self.epsg = epsg

    def assign_geometry(self):

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
        srs.ImportFromEPSG(self.epsg)
        outdataset.SetProjection(srs.ExportToWkt())
        count=0
        for arr in img_array:
            count+=1
            outband=outdataset.GetRasterBand(count)
            outband.WriteArray(arr)
            outband=None
        
        return(outfile)