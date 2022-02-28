# The Automated Polynya Identification Tool (APIT).

## Quick script to image reference:
|Step | Script | Input |
|----|---------|:-------:|
|1 Acquisition |           |       |
|1.1    | 01_extract_urls.py            |       |
|1.2    | 02_LPDAAC_download.py         |       |
|1.3    | 03_rename_files.py            |       |
|2 Pre-process |           |       |               |
|2.1    | 01_modis_preprocess.py        | 01_file.jpg/.hdf |
|2.2    | 02_modis_data_fusion.py       | 02_file.tif |
|3      |       |       |               |
|3.1    | 01_threshold_classification.py| 02_file.tif |
|3.2    | 02_noise_removal.py           | 03a_file.tif |
|3.3    | 03_filter.py                  | 03c_file.tif |
|4      |       |       |
|4.1    | 01_netcdf.py      | 03f_file.tif      |
