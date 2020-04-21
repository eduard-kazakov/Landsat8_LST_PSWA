# Landsat8_LST_PSWA
Open and fully automated python realization of Practical Split-Window Algorithm (PSWA) for Land Surface Temperature (LST) retrieval from Landsat-8/TIRS data. Covariance-variance ratio (SWCVR) method is used for water vapour content calculations. FMASK is used for cloud/water masks obtaining.

**Important note**: python-fmask (http://www.pythonfmask.org/en/latest/) module must be installed and available with command-line interface

## References

Theoretical basis for implemented methods are available in following papers:

* Du, Chen; Ren, Huazhong; Qin, Qiming; Meng, Jinjie; Zhao, Shaohua. 2015. "A Practical Split-Window Algorithm for Estimating Land Surface Temperature from Landsat 8 Data." Remote Sens. 7, no. 1: 647-665
* Huazhong Ren, Chen Du, Qiming Qin, Rongyuan Liu, Jinjie Meng, and Jing Li. "Atmospheric Water Vapor Retrieval from Landsat 8 and Its Validation." 3045--3048. IEEE, 2014.
* Ren, H., Du, C., Liu, R., Qin, Q., Yan, G., Li, Z. L., & Meng, J. (2015). Atmospheric water vapor retrieval from Landsat 8 thermal infrared images. Journal of Geophysical Research: Atmospheres, 120(5), 1723-1738

## Requirements

gdal, numpy, python-fmask (must be available in command-line interface)

Optional requirement: SREMPy-landsat (https://github.com/eduard-kazakov/SREMPy-landsat). If installed, additional option of LSE calculation is available.

## Recommended environment preparation

I recommend to use conda (miniconda build) (https://docs.conda.io/en/latest/miniconda.html):

```bash
conda create --name lst python=3.8.2
conda install -c conda-forge -n lst python-fmask
conda activate lst
pip install git+https://github.com/eduard-kazakov/SREMPy-landsat
```

## Installation

```bash
pip install git+https://github.com/eduard-kazakov/Landsat8_LST_PSWA
```

## For what this algorithm is developed?

Land surface temperature is very important geophysical parameter. Many Earth Observation Satellites are equiped with thermal sensors, which are able to measure surface temperature. But there are a big problem: with remote sensing we can directly measure only **brightness temperature** (or radiative temperature). This parameter is not equal to **physical surface temperature**, and not taking into account surface material, atmospherical conditions and so on (and this is VERY important when you compare temperatures or calculate some temperature-based products). There are many approaches to retrieve real surface LST from measured on satellite brightness temperature, but most of them require a lot of additional information about atmosphere. Nice review of such methods could be found here:
* Jiménez-Muñoz J. C. et al. Land surface temperature retrieval methods from Landsat-8 thermal infrared sensor data //IEEE Geoscience and remote sensing letters. – 2014. – Т. 11. – №. 10. – С. 1840-1843.
* Yu X., Guo X., Wu Z. Land surface temperature retrieval from Landsat 8 TIRS—Comparison between radiative transfer equation-based method, split window algorithm and single channel method //Remote sensing. – 2014. – Т. 6. – №. 10. – С. 9829-9852.

For operational purposes we usually want to have fully automated algorithms, which are not have external data dependencies. Practical Split-Window algorithm (Du, Qiu, Meng, Zhao, 2015) was developed for such cases and have satisfactory quality. Water vapour content (required for Split-Window algorithms) is obtained with automated SWCVR method, based on values ratio in 10 and 11 bands. This python implementation allows to calculate physicall LST without any external data and in fully automated mode. Python-fmask algorithm is used for masking clouds and water.      

## How to use

All options are set with class object (LSTRetriever) initialization. Initialization arguments are:
* **metadata_file** - full path to Landsat 8 (L1C) dataset metadata file (*MTL.txt)
* **LSE_mode** - is option how to calculate Land Surface Emissivity (LSE), very important parameter for LST. Supported variants: 
    * auto-ndvi-raw : LSE is calculated from NDVI reflectace values (without atmosperic correction)
    * auto-ndvi-srem : LSE is calculated from NDVI reflectace values (with SREM atmosperic correction). SREMPy-landsat module must be installed for this option.
    * from-glc : LSE is obtained from land cover maps (GLC) with look-up table. FROM-GLC project is used as land cover source (http://data.ess.tsinghua.edu.cn/)
    * external : LSE is obtained from external file (if you can calculate LSE elsewhere)
* **LSE_file** - is used only if LSE is **from-glc** or **external**
    * if **from-glc** mode is set, this file must be GeoTiff with land cover map (with original FROM-GLC values)
    * if **external** mode is set, this file must be 2-band GeoTiff with LSE for B10 in first band and LSE for B11 in second band
* **temp_dir** - temporary directory path.
* **window_size** - split window size. Default value is 15 (11x11 pixels). You can experiment with huge values for better runtime performance. Small window size leads to long processing time.
* **angles_file** - is used only if LSE is **auto-ndvi-srem**. This is path to *ANG.txt file from standart Landsat 8 L1C dataset
* **usgs_utils** - is used only if LSE is **auto-ndvi-srem**. This is path to compiled executable of l8_angles util by USGS (https://www.usgs.gov/land-resources/nli/landsat/solar-illumination-and-sensor-viewing-angle-coefficient-files). Further information about SREMPy-landsat is available here: https://github.com/eduard-kazakov/SREMPy-landsat.
* **cygwin_bash_exe_path** - only for windows users who wants to use SREM for LSE.

Simplest option is to use **auto-ndvi-raw**, but it also worst in quality.

After initialization you can directly retrieve products with following methods:
* **get_lst_as_gtiff (output_path)** - main function, doing all work and generating output LST GeoTiff
* **get_lst_array** - same, but without GeoTiff creating. Returns array, so you can use it for next processing steps
* **get_cwv_as_array** - get array with water vapour content (in same extent as original L8 scene)
* **get_b10_b11_brightness_temp_arrays** - get two arrays with brightness temperatures for B10 and B11 (in same extent as original L8 scene)
* **__save_array_to_gtiff (output_path)** - you can use this method to save any array (generated with described methods) to GeoTiff.

You can also modify class properties, affecting cloud masking (after initialization, before running lst calculation):
* cloudbufferdistance (default 30)
* cloudprobthreshold (default 75.0)
* shadowbufferdistance (default 0)

## Examples

1. With **auto-ndvi-raw** mode:

```python         
from L8LST_PSWA.LSTRetriever import LSTRetriever

lst_retriever = LSTRetriever(metadata_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_MTL.txt',
                            LSE_mode='auto-ndvi-raw',
                            temp_dir='/.../temp_dir',
                            window_size=15)

lst_retriever.cloudprobthreshold = 60
lst_retriever.get_lst_as_gtiff('/.../LST_20180512.tif')
```

2. With **from-glc** mode:

```python         
from L8LST_PSWA.LSTRetriever import LSTRetriever

lst_retriever = LSTRetriever(metadata_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_MTL.txt',
                            LSE_mode='from-glc',
                            LSE_file='/.../Data/glc_nw.tif',                            temp_dir='/.../temp_dir',
                            window_size=15)

lst_retriever.get_lst_as_gtiff('/.../LST_20180512.tif')
```


3. With **auto-ndvi-srem** mode:

```python         
from L8LST_PSWA.LSTRetriever import LSTRetriever

lst_retriever = LSTRetriever(metadata_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_MTL.txt',
                            LSE_mode='auto-ndvi-srem',
                            angles_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_ANG.txt',
                            usgs_utils='/.../L8_ANGLES_2_7_0/l8_angles/l8_angles',
                            temp_dir='/.../temp_dir',
                            window_size=15)

lst_retriever.get_lst_as_gtiff('/.../LST_20180512.tif')
``` 

4. With **auto-ndvi-srem** mode (windows users):

```python         
from L8LST_PSWA.LSTRetriever import LSTRetriever

lst_retriever = LSTRetriever(metadata_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_MTL.txt',
                            LSE_mode='auto-ndvi-srem',
                            angles_file='/.../Data/LC08_L1TP_185018_20180512_20180517_01_T1/LC08_L1TP_185018_20180512_20180517_01_T1_ANG.txt',
                            usgs_utils='/.../L8_ANGLES_2_7_0/l8_angles/l8_angles.exe',
                            temp_dir='/.../temp_dir',
                            window_size=15,
                            cygwin_bash_exe_path='C:/cygwin64/bin/bash.exe')

lst_retriever.get_lst_as_gtiff('/.../LST_20180512.tif')
``` 

## Contacts

Feel free to contact me: ee.kazakov@gmail.com
