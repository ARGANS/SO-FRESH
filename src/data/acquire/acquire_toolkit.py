#!/usr/bin/env/python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing data acquisition for the Automated Polynya Identification Tool.
"""

# Package loader # 
import os, sys
from bs4 import BeautifulSoup
from getpass import getpass
from netrc import netrc
import requests
from tqdm import tqdm

# Toolkit loader # 
from ..adjust import adjust_toolkit as adj_tk

class acquire_modis:
    def __init__(self, data_directory, acquire_dates , product, aoi):
        self.data_directory = data_directory
        self.acquire_dates = acquire_dates
        self.product = product
        self.aoi = aoi
        

    def product_url(self):

        """ Identifies product requested for download and specifies the URL. """
        print(f"Pulling {self.product} product link...\n")
        if self.product == "MYD09GA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"
        elif self.product == "MYDTBGA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/"

        return(url)

    def url_extractor(self, url):

        """ Pulls the URL of all products which have fit the criteria and creates their out directory """

        hmin, hmax, vmin, vmax = int(self.aoi[0]), int(self.aoi[1]), int(self.aoi[2]), int(self.aoi[3])
        if str(self.aoi) == "(14, 24, 15, 17)":roi = "antarctica"
        elif str(self.aoi) == "(13, 23, 0, 2)":roi = "arctic"
        item_links, outdirs=[],[]
        print("Pulling links of all products that fit input parameters:")
        for date in tqdm(self.acquire_dates):
            dpage = requests.get(url+date).text
            soup = BeautifulSoup(dpage, "lxml")
            if self.product == "MYD09GA":
                for link in soup.select("a[href$='.jpg']"):
                    jpgPath = (url+date+link.get("href"))
                    h = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][1:3]
                    v = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][4:]
                    if int(h) >= hmin and int(h)<= hmax:
                        if int(v) >= vmin and int(v) <= vmax:
                            outdir = self.data_directory+"02_data/MODIS/"+self.product+"_061/01_tiles/"+"h"+h+"v"+v+"/"+str(date).replace(".", "/")
                            if not os.path.isdir(outdir): os.makedirs(outdir)
                            prod_lnk = ("%s\n"%jpgPath)
                            outdirs.append(outdir)
                            item_links.append(prod_lnk)

            if self.product == "MYDTBGA":
                ################# CHECK IF THIS IS REAL - IF TRUE DON'T DOWNLOAD #####################
                print(self.data_directory+"02_data/MODIS/"+self.product+"_006/02_mosaic/"+str(date).replace(".", "/")+"02_"+self.product+"_006_"+str(date).replace(".", "")[:-1]+"_"+roi+".tif")
                print(str(self.aoi))
                ### also sort adjust process for this data ###
                #sys.exit()


                sys.exit()
                for link in soup.select("a[href$='.hdf']"):
                    hdfPath = (url+date+link.get("href"))
                    h = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][1:3]
                    v = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][4:]
                    if int(h) >= hmin and int(h)<= hmax:
                        if int(v) >= vmin and int(v) <= vmax:
                            outdir = self.data_directory+"02_data/MODIS/"+self.product+"_006/01_tiles/"+"h"+h+"v"+v+"/"+str(date).replace(".", "/")
                            if not os.path.isdir(outdir): os.makedirs(outdir)
                            prod_lnk = ("%s\n"%hdfPath)
                            outdirs.append(outdir)
                            item_links.append(prod_lnk)

        return(item_links, outdirs)

class lpdaac:
    @staticmethod
    def earthdata_authentication():

        """ Login for Earth data. """

        urs="urs.earthdata.nasa.gov"
        prompts = ['Enter NASA Earthdata Login Username \n(or create an account at urs.earthdata.nasa.gov): ','Enter NASA Earthdata Login Password: ']
        try: # Determines if netrc file exists with Earth Data credentials.
            netrcdir=os.path.expanduser("~/.netrc")
            netrc(netrcdir).authenticators(urs)[0]

        except FileNotFoundError: # Creates netrc and prompts user for Earthdata login details
            homedir=os.path.expanduser("~")
            Popen('touch {0}.netrc | chmod og-rw {0}.netrc | echo machine {1} >> {0}.netrc'.format(homedir + os.sep, urs), shell=True)
            Popen('echo login {} >> {}.netrc'.format(getpass(prompt=prompts[0]), homedir + os.sep), shell=True)
            Popen('echo password {} >> {}.netrc'.format(getpass(prompt=prompts[1]), homedir + os.sep), shell=True)

        except TypeError: # Edit netrc file is it is not setup for Earthdata.
            homedir = os.path.expanduser("~")
            Popen('echo machine {1} >> {0}.netrc'.format(homedir + os.sep, urs), shell=True)
            Popen('echo login {} >> {}.netrc'.format(getpass(prompt=prompts[0]), homedir + os.sep), shell=True)
            Popen('echo password {} >> {}.netrc'.format(getpass(prompt=prompts[1]), homedir + os.sep), shell=True)
        # Gives the user one minute to submit username and password.
        tries = 0
        while tries < 30:
            try:
                netrc(netrcdir).authenticators(urs)[2]
            except:
                time.sleep(2.0)
            tries += 1

        return(netrcdir, urs)

    @staticmethod
    def workspace_setup(data_directory, url):

        """ Preparation for data download. """
        
        if isinstance(url, str):
            file_lst = [url]
        if data_directory[-1] != "/" and data_directory[-1] != "\\":
            data_directory=data_directory.strip("'").strip('"')+os.sep
    
    @staticmethod
    def execute_download(data_directory, links, dirs, resample=None, adjust=None):
        
        """ Execute data download from Earth data. """

        imgs_4_mosaic = []
        netrcdir, urs = lpdaac.earthdata_authentication()
        if adjust == True: print("Downloading and adjusting...\n")
        else: print("Downloading...\n")
        for l, d in tqdm(zip(links, dirs), total=len(links)):
            lpdaac.workspace_setup(data_directory, l)
            outfile=d+"01_"+str(l.split("/")[-1].strip())
            if os.path.exists(outfile) or os.path.exists((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif"):
                imgs_4_mosaic.append((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif")
                pass   
            else:
                with requests.get(l.strip(), verify=False, stream=True, auth=(netrc(netrcdir).authenticators(urs)[0], netrc(netrcdir).authenticators(urs)[2])) as response:
                    if response.status_code != 200:
                        print("{} not downloaded. Verify that your username and password are correct in {}".format(l.split('/')[-1].strip(), netrcdir))
                    else:
                        response.raw.decode_content = True
                        content = response.raw
                        with open(outfile, 'wb') as d:
                            while True:
                                chunk = content.read(16 * 1024)
                                if not chunk:
                                    break
                                d.write(chunk)
                        print(f"Downloaded file: {l}")
            if adjust == True:
                if os.path.exists((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif"):
                    imgs_4_mosaic.append((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif")
                    continue
                elif os.path.exists(outfile) and not os.path.exists((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif"):
                    if "MYD09GA" in outfile:
                        OPT_adj = adj_tk.myd09ga_adjust(outfile)
                        rpjct = OPT_adj.reproject()
                        if not resample == "False":
                            rename = os.path.split(rpjct)[0]+"/"+os.path.split(rpjct)[1][:2]+"a"+os.path.split(rpjct)[1][2:]
                            os.rename(rpjct, rename)
                            xres, yres = resample[0], resample[1]
                            rsampled = adj_tk.myd_adjust.resample(rename, xres, yres)
                            os.remove(outfile)
                            os.remove(rename)
                            imgs_4_mosaic.append(rsampled)
                        elif resample == "False":
                            print("Appending reprojected")
                            imgs_4_mosaic.append(rpjct)
                        else:pass
                    elif "MYDTBGA" in outfile:
                        print("MYDTBGA processing to be set-up.")
                        # extract, normalise, assign geometry and resample
                        sys.exit()

            else:continue
        return(imgs_4_mosaic)