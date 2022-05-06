#!/usr/bin/env/python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing data acquisition for the Automated Polynya Identification Tool.
"""
# modules
import os, sys
from bs4 import BeautifulSoup
from datetime import date, timedelta
from getpass import getpass
import glob
import itertools
from netrc import netrc
import pandas as pd
import requests
import shutil
import time
from tqdm import tqdm
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+"/adjust")
import toolkit_ as tk_


class lpdaac_download():
    def __init__(self, data_folder, sdate, edate, product, hmin, hmax, vmin, vmax):
        self.data_folder=data_folder
        self.sdate=sdate
        self.edate=edate
        self.product=product
        self.hmin=hmin
        self.hmax=hmax
        self.vmin=vmin
        self.vmax=vmax

    def dates_builder(self):

        """ Identify all dates between start and end entries. """

        dates=[]

        syear=int(self.sdate[0:4])
        smonth=int(self.sdate[5:7])
        sday=int(self.sdate[8:10])

        d0 = date(syear, smonth, sday)

        eyear=int(self.edate[0:4])
        emonth=int(self.edate[5:7])
        eday=int(self.edate[8:10])

        d1 = date(eyear, emonth, eday)
        delta=d1-d0
        for i in range(delta.days+1):
            day=d0+timedelta(days=i)
            dates.append(str(day).replace("-", ".")+"/")

        return(dates)

    def product_selector(self):

        """ Identifies product requested for download and specifies the URL. """

        if self.product == "MYD09GA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"
        elif self.product == "MYDTBGA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/"

        return(url)

    def url_extractor(self):

        """ Pulls the URL of all products which have fit the criteria and creates their out directory """

        url=self.product_selector()
        prods_links, outdirs=[],[]
        for date in tqdm(self.dates_builder()):
            dpage = requests.get(url+date).text
            soup = BeautifulSoup(dpage, "lxml")
            if self.product == "MYD09GA":
                for link in soup.select("a[href$='.jpg']"):
                    jpgPath = (url+date+link.get("href"))
                    h = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][1:3]
                    v = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][4:]
                    if h >= self.hmin and h<= self.hmax:
                        if v >= self.vmin and v <= self.vmax:
                            prod_lnk = ("%s\n"%jpgPath)
                            # execute download here
                            ##### CODE TO BE ADDED FOR OPTICAL DOWNLOAD AND PROCESSING ####
            if self.product == "MYDTBGA":
                for link in soup.select("a[href$='.hdf']"):
                    hdfPath = (url+date+link.get("href"))
                    h = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][1:3]
                    v = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][4:]
                    if h >= self.hmin and h<= self.hmax:
                        if v >= self.vmin and v <= self.vmax:
                            outdir = self.data_folder+"MODIS/"+self.product+"_006/01_tiles/"+"h"+h+"v"+v+"/"+str(date).replace(".", "/")
                            if not os.path.isdir(outdir): os.makedirs(outdir)
                            prod_lnk = ("%s\n"%hdfPath)
                            outdirs.append(outdir)
                            prods_links.append(prod_lnk)

        return(prods_links, outdirs)


    def workspace_setup(self, url):

        """ Preparation for data download. """
        
        if isinstance(url, str):
            file_lst = [url]
        if self.data_folder[-1] != "/" and self.data_folder[-1] != "\\":
            self.data_folder=self.data_folder.strip("'").strip('"')+os.sep

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


    def download(self, links, dirs):

        """ Execute data download from Earth data. """

        netrcdir, urs = self.earthdata_authentication()
        for l, d in tqdm(zip(links, dirs), total=len(links)):
            self.workspace_setup(l)
            outfile=d+"01_"+str(l.split("/")[-1].strip())
            if os.path.exists(outfile) or os.path.exists((os.path.split(outfile))[0]+"/02_"+(os.path.split(outfile)[1])[3:-3]+"tif"):
                continue
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
                
                extract=tk_.MYDTBGA_preprocess(outfile).extract_MYDTBGA()
                normalise=tk_.MYDTBGA_preprocess((extract+"_4.tif")).normalise()
                os.remove(outfile)
                shutil.rmtree(os.path.split(extract)[0], ignore_errors=True)
        
        return(normalise)

class scan_database():
    def __init__(self, data_folder, sdate, edate, product, hmin, hmax, vmin, vmax):
        self.data_folder=data_folder
        self.sdate=sdate
        self.edate=edate
        self.product=product
        self.hmin=hmin
        self.hmax=hmax
        self.vmin=vmin
        self.vmax=vmax

    def dates_builder(self):

        """ Identify all dates between start and end entries. """

        dates=[]

        syear=int(self.sdate[0:4])
        smonth=int(self.sdate[5:7])
        sday=int(self.sdate[8:10])

        d0 = date(syear, smonth, sday)

        eyear=int(self.edate[0:4])
        emonth=int(self.edate[5:7])
        eday=int(self.edate[8:10])

        d1 = date(eyear, emonth, eday)
        delta=d1-d0
        for i in range(delta.days+1):
            day=d0+timedelta(days=i)
            dates.append(str(day).replace("-", "/")+"/")

        return(dates)

    def tiles_builder(self, csv=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../", "aux_files/modis_sinusoidal_tiles.csv"))):
        
        """ Identify all tile posibilities from entries. """
        
        h = [self.hmin+i for i in range((self.hmax - self.hmin)+1)]
        v = [self.vmin+i for i in range((self.vmax - self.vmin)+1)]
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

    def file_iterator(self, dates=None,tiles=None,root=None, version=None):

        """ Itterates through specified range of dates and tiles, checks present files and downloads or processes, based on what is available. """

        for t in tiles:
            h, v = t[1:-3], t[4:]
            for d in dates:
                fp = root+t+"/"+d
                if not os.path.isdir(fp): os.makedirs(fp)
                if os.path.isdir(fp):
                    files=glob.glob(fp+"*")
                    if len(files) == 0:
                        print(f"Searching for files for h{h}v{v} on {d}...")
                        links, outdirs = self.individual_url_extractor(ind_dates=d, h=h, v=v, version=version)
                        if links == None and outdirs == None:
                            dwnld_output=None
                            print("No files found, passing...")
                        else:    
                            print("Downloading and processing...\n")
                            dwnld_output = lpdaac_download(self.data_folder, self.sdate, self.edate, self.product, self.hmin, self.hmax, self.vmin, self.vmax).download([links], [outdirs])
                    elif (all(os.path.basename(item)[0:2] == "01" for item in files)) and len(files) == 1:
                        # If there is just the raw files - delete & re-download.
                        # The extraction & normalisation has failed in prior download step.
                        print(f"Found an issue with tile h{h}v{v} on {d}, re-downloading and processing...")
                        os.remove(files[0])
                        links, outdirs = self.individual_url_extractor(ind_dates=d, h=h, v=v, version=version)
                        dwnld_output = lpdaac_download(self.data_folder, self.sdate, self.edate, self.product, self.hmin, self.hmax, self.vmin, self.vmax).download([links], [outdirs])
                    elif (all(os.path.basename(item)[0:2] == "01" for item in files)) and len(files) == 2:
                        # If there is the raw files & extraction directory (ED), but ED is empty, the raw file is corrupt.
                        # Therefore, delete the raw file and ED, and re-download and process. 
                        f = [f for f in files if os.path.basename(f).startswith("01a") == True][0]
                        if len(os.listdir(f)) < 7:
                            print(f"Found a corrupt file for for tile h{h}v{v} on {d}, re-downloading and processing...")
                            shutil.rmtree(f)
                            os.remove(files[0])
                            links, outdirs = self.individual_url_extractor(ind_dates=d, h=h, v=v, version=version)
                            dwnld_output = lpdaac_download(self.data_folder, self.sdate, self.edate, self.product, self.hmin, self.hmax, self.vmin, self.vmax).download([links], [outdirs])
                        elif len(os.listdir(f)) == 7 and any(fname.endswith('_4.tif') for fname in os.listdir(f)):
                            # If the extraction directory happens to have the files in, then normalise. 
                            print(f"Found un-processed files for tile h{h}v{v} on {d}, normalising...")
                            tk_.MYDTBGA_preprocess(f+"/"+[b for b in os.listdir(f) if b.endswith("_4.tif") == True][0]).normalise()
                            # Remove raw file and extracted hdf files.
                            os.remove([f for f in files if os.path.basename(f).startswith("01_") == True][0])
                            shutil.rmtree(f)
                    elif (all(os.path.basename(item)[0:2] == "02" for item in files)):
                        # If required file exists, then continue.
                        continue
        return(dwnld_output)

              
    def individual_url_extractor(self, ind_dates=None, h=None, v=None, version=None):

        """ Pulls the URL of a specified product which fits the criteria and creates its output directory """

        url= self.product_selector()
        prods_links, outdirs=[],[]
        for date in ind_dates:
            date = (str(ind_dates).replace("/", ".")[:-1]+"/")
            dpage = requests.get(url+date).text
            soup = BeautifulSoup(dpage, "lxml")

            if self.product == "MYD09GA":
                for link in soup.select("a[href$='.jpg']"):
                    print("here")
                    jpgPath = (url+date+link.get("href"))
                    h_Path = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][1:3]
                    v_Path = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][4:]
                    if h_Path != h and v_Path != v:
                        continue
                    elif h_Path == h and v_Path == v:
                        outdir = self.data_folder+"MODIS/"+self.product+version+"/01_tiles/"+"h"+h_Path+"v"+v_Path+"/"+str(date).replace(".", "/")
                        if not os.path.isdir(outdir): os.makedirs(outdir)
                        prod_lnk = ("%s\n"%jpgPath)
                        return(prod_lnk, outdir)
                return(None, None)

            if self.product == "MYDTBGA":
                for link in soup.select("a[href$='.hdf']"):
                    hdfPath = (url+date+link.get("href"))
                    h_Path = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][1:3]
                    v_Path = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][4:]
                    if h_Path != h and v_Path != v:
                        continue
                    elif h_Path == h and v_Path == v:
                        outdir = self.data_folder+"MODIS/"+self.product+"_006/01_tiles/"+"h"+h_Path+"v"+v_Path+"/"+str(date).replace(".", "/")
                        if not os.path.isdir(outdir): os.makedirs(outdir)
                        prod_lnk = ("%s\n"%hdfPath)
                        return(prod_lnk, outdir)
                return(None, None)

    def product_selector(self):

        """ Identifies product requested for download and specifies the URL. """

        if self.product == "MYD09GA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"
        elif self.product == "MYDTBGA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/"

        return(url)
