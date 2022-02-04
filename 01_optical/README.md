This is the development stage of the Automated Polynya Identification Tool (APIT) v1, using MODIS (MYD09GA) data. Information surrounding this tool will be documented on the [wiki](https://argans.atlassian.net/wiki/spaces/SO/overview?homepageId=995393752 "SO-Fresh Wiki").

## Table of contents:
* [01_download](https://github.com/ARGANS/SO-FRESH/tree/development#01_download)
* [02_preprocess](https://github.com/ARGANS/SO-FRESH/tree/development#02_preprocess)
* [03_classification](https://github.com/ARGANS/SO-FRESH/tree/development#03_classification)
* [04_identify](https://github.com/ARGANS/SO-FRESH/tree/development#04_identify)
* [05_filter](https://github.com/ARGANS/SO-FRESH/tree/development#05_filter)
* [06_compare](https://github.com/ARGANS/SO-FRESH/tree/development#06_compare)

#### Quick-reference script-->input image.
|#    | Script        | Quick-input reference     |
| --- | ------------- |:-------------:|
| 1  | 02_preprocess/modis_prjct_jpg.py     |       01*.jpg     |
| 2  | 03_classification/threshold_classifier.py       |       02*.tif      |
| 3  | 04_identify/polynya_identification_heatmap.py       |       03*.tif      |
| 4  | 05_filter/measure_mask_area_v2.py       |       05*.tif      |
| 5  | 06_compare/assess_over_time.py       |       08*.tif      |

***

# 01_download
MODIS imagery is downloaded from the [The Land Processes Distributed Active Archive Centre](https://lpdaac.usgs.gov/ "LPDAAC"), where "True-Colour Images" of NASA products are able to be downloaded. 

The product used for this tool is MODIS MYD09GA, where images of each MODIS tile:
* Contain RGB band information.
* Are ~ 30 - 40 kb in size.
* Rich archive dating from 2002 - present. 
* Available in version [006](https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/ "MYD09GA.006") and [061](https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/ "MYD09GA.061").

## 1.1 extract_urls.py
For product download, exact urls to individual products are required. This script extracts the URLS to each product through webscraping and puts them in a text file.

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

## 1.2 DAAC_data_download.py

Using the text file created prior, uses all available URLS to download those products to the desired directory.

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

## 1.3 Rename_Files.py

Adds '01' to the begining of the image filename.

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| inDir     |       -d      | Filepath to the directory with all the files which are to be renamed. |

###### Example:
```
python SO-FRESH/01_download/Rename_Files.py -inDir download_imagery/tile/date/file.jpg
```
***

## 02_preprocess
Imagery downloaded from the LPDAAC website comes unprojected and in a JPEG format. In this directory is a CSV which contains the co-ordinates of all MODIS sinusoidal tiles.

## 2.1 modis_projct_jpg.py

Based on tile reference (taken from image name), project image using co-ordinates in CSV. The output is the same directory as the input, with a changed name.

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| input-img     |       -i      | Input image to be projected. |

```
python SO-FRESH/02_preprocess/modis_prjct_jpg.py -i download_imagery/tile/date/01_file.jpg
```
***

## 03_classification
Undertake a threshold classification based on a set value (v1).

## 3.1 threshold_classifier.py

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| input-img     |       -i      | Input image to be classified. |
| threshold     |       -t      | Threshold value to use. |

```
python SO-FRESH/03_classification/threshold_classifier.py -i download_imagery/tile/date/02_file.tif -t 180
```
***

## 04_identify
**This section requires a rename**
From the classification, two stages are produced in this section:
* Heatmap is produced from a moving window
* Followed by a noise removal stage, removing any small data surrounding the open-water areas.

## 4.1 polynya_identification_heatmap.py

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| input-img     |       -i      | Input classification image to be filtered. |

```
python SO-FRESH/04_identify/polynya_identification_heatmap.py -i download_imagery/tile/date/03_file_classification.tif 
```
***

## 05_filter
From the classification, potential areas are vectorized to calculate area and distance from land. Based on this criteria, data is masked and removed.

## measure_mask_area_v2.py

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| input-img     |       -i      | Input noise removal'd image. |

```
python SO-FRESH/05_filter/measure_mask_area_v2.py -i download_imagery/tile/date/05_file_classification_noise_removal.tif 
```
***
## 06_compare
Assess the occurance of open-water formations over a given period of time. 

## assess_over_time.py

| Inputs        | Shorthand     | What is it?  |
| ------------- |:-------------:| ------------:|
| time-start     |       -s      | Date to start from (YYYY/MM/DD). |
| tile     |       -t      | Filepath to tile directory of interest. |
| days     |       -d      | Number of days to look at from the start-date. |
| time-end     |       -e      | Look up to this period. |
| save     |       -save      | Save image in current directory as 'h##v##_YYYYMMDD_YYYYMMDD.tif'. |
| output-img     |       -o      | Output image is specified directory and filename. |

###### Extra information:
* Start date and tile are **required**, alongside either number of days or time-end.
***