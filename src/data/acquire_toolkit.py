#!/usr/bin/env/python3
# @James Hickson | Argans UK | jhickson@argans.co.uk

import os, sys
from bs4 import BeautifulSoup
from datetime import date, timedelta
from getpass import getpass
from netrc import netrc
import requests
import shutil
import time
from tqdm import tqdm
import toolkit as tk

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
        if self.product == "MYD09GA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"
        elif self.product == "MYDTBGA":
            response = requests.get("https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/").text
            url = "https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/"

        return(url)

    def url_extractor(self):
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
        if isinstance(url, str):
            file_lst = [url]
        if self.data_folder[-1] != "/" and self.data_folder[-1] != "\\":
            self.data_folder=self.data_folder.strip("'").strip('"')+os.sep

    @staticmethod
    def earthdata_authentication():
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
        netrcdir, urs = self.earthdata_authentication()
        for l, d in zip(links, dirs):
            self.workspace_setup(l)
            outfile=d+"01_"+str(l.split("/")[-1].strip())
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
            extract=tk.modis_preprocess(outfile).extract_MYDTBGA()
            normalise=tk.modis_preprocess((extract+"_4.tif")).normalise()
            os.remove(outfile)
            shutil.rmtree(os.path.split(extract)[0], ignore_errors=True)
            sys.exit()
