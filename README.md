This is the development stage of the Automated Polynya Identification Tool (APIT) v1, using MODIS (MYD09GA) data. Information surrounding this tool will be documented on the [wiki](https://argans.atlassian.net/wiki/spaces/SO/overview?homepageId=995393752 "SO-Fresh Wiki").

***
<details>
    <summary>01_download</summary>

    MODIS imagery is downloaded from the [The Land Processes Distributed Active Archive Centre](https://lpdaac.usgs.gov/ "LPDAAC"), where "True-Colour Images" of NASA products are able to be downloaded. 

    The product used for this tool is MODIS MYD09GA, where images of each MODIS tile:
    * Contain RGB band information.
    * Are ~ 30 - 40 kb in size.
    * Rich archive dating from 2002 - present. 
    * Available in version [006](https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/ "MYD09GA.006") and [061](https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/ "MYD09GA.061").
<details>
    <summary>1.1 extract_urls.py</summary>
    
    | Inputs        | Shorthand     | What is it?  |
    | ------------- |:-------------:| ------------:|
    | startDate     |       -s      | Start date (YYYY-MM-DD). |
    | endDate       |       -e      | End date (YYYY-MM-DD). |
    | outpath       |       -o      |    Path to where the textfile is saved. |
    | version       |       -v      |    MODIS data version (006 or 061). |

    ###### Example:
    ```
    python SO-FRESH/01_download/extract_urls.py -s 2017-01-01 -e 2017-12-31 -o download_text/ -v 006
    ```
    ###### Extra information:
    * Line 97 - The final value is the tile of interest - this should be modified based on v tile of interest.
<details>

<details>
    <summary>1.2 DAAC_data_download.py</summary>

    | Inputs        | Shorthand     | What is it?  |
    | ------------- |:-------------:| ------------:|
    | directory     |       -dir      | Specification of output directory. |
    | files       |       -f      | Filepath to textfile created in previous step. |

    ###### Example:
    ```
    python SO-FRESH/01_download/DAAC_data_download.py -dir download_imagery/ -f download_text/imagery.txt
    ```
    ###### Extra information:
    * Source: https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_python/browse
    * Login details are required in the '.netrc' file in the following format:
    ```
    machine urs.earthdata.nasa.gov
    login jhickson
    password password123
    ```
<details>
<details>

***

## 02_preprocess

## 03_classification

## 04_identify

## 05_filter

## 06_compare