# The Automated Polynya Identification Tool (APIT).

## Quick script to image reference:
|Step | Script | Input |
|----|---------|:-------:|
|1.1    | 01_extract_urls.py            | ---   |
|1.2    | 02_LPDAAC_download.py         | ---   |
|1.3    | 03_rename_files.py            | ---   |
|2.1    | 01_modis_preprocess.py        | 01_file.jpg/.hdf |
|2.2    | 02_modis_data_fusion.py       | 02_file.tif |
|3.1    | 01_threshold_classification.py| 02_file.tif |
|3.2    | 02_noise_removal.py           | 03a_file.tif |
|3.3    | 03_filter.py                  | 03c_file.tif |
|4.1    | 01_netcdf.py                  | 03f_file.tif      |
**
# 1. Data acquisition and download.
MODIS data is acquired from [The Land Processes Distributed Active Archive Centre](https://lpdaac.usgs.gov/ "LPDAAC").
From this archive centre, current data available for download is:
* MYD09GA
* MYDTBGA

MYD09GA = optical imagery (RGB)
MYDTBGA = thermal imagery (band 32)

## 1.1 01_extract_urls.py
Extract full links for product download, through the use of web-scraping, delivering a textfile output which each link as a new line.
For this, it is recommended that you do year time-stamp intervals (this makes the next step easier). 

| Inputs        | Shorthand     | What is it?  | Notes |
| ------------- |:-------------:| ------------:|
| time-start     |       -s      | Start date (YYYY-MM-DD). |
| time-end       |       -e      | End date (YYYY-MM-DD). |
| outpath       |       -o      |    Path to where the textfile is saved. |
| version       |       -v      |    006 or 061 | V006 to be discontinued |
| data-type |           -t |        Optical or thermal
|horizontal-minimum|    -hmin | Horizontal minimum MODIS sinusoidal tile|
|horizontal-maximum|    -hmax | Horizontal maximum MODIS sinusoidal tile|
|vertical-minimum|      -vmin | Vertical minimum MODIS sinusoidal tile |
|vertical-maximum|      -vmax | Vertical maximum MODIS sinusoidal tile |

## 1.2 02_LPDAAC_download.py
Initiate data download. If the date range is too broad, the download may not be able to begin due to over-requesting from the archive centre. Shell scripting this is a recommendation.

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| directory     |       -dir      | Specification of output directory. |
| files       |       -f      | Filepath to textfile created in previous step. |

Ensure your .netrc file in your home directory contains the following:
```
machine urs.earthdata.nasa.gov
login jhickson
password password123
```


