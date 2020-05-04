# -*- coding: utf-8 -*-

# Practical split-window algorithm + covariance-variance ratio (SWCVR) method for Landsat 8 imagery
#
# References:
# Du, Chen; Ren, Huazhong; Qin, Qiming; Meng, Jinjie; Zhao, Shaohua. 2015. "A Practical Split-Window Algorithm for Estimating Land Surface Temperature from Landsat 8 Data." Remote Sens. 7, no. 1: 647-665
# Huazhong Ren, Chen Du, Qiming Qin, Rongyuan Liu, Jinjie Meng, and Jing Li. "Atmospheric Water Vapor Retrieval from Landsat 8 and Its Validation." 3045--3048. IEEE, 2014.
# Ren, H., Du, C., Liu, R., Qin, Q., Yan, G., Li, Z. L., & Meng, J. (2015). Atmospheric water vapor retrieval from Landsat 8 thermal infrared images. Journal of Geophysical Research: Atmospheres, 120(5), 1723-1738
#
# Eduard Kazakov, 2020 | ee.kazakov@gmail.com

import os
import math

import gdal
import numpy as np

FMASK_EXECUTABLE_PATH = 'fmask_usgsLandsatStacked.py'

class CalibrateLandsatBand():
    def __init__(self, band_file, metadata_file):
        self.metadata_reader = LandsatMetadataReader(metadata_file)

        self.metadata = self.metadata_reader.metadata

        self.band_metadata = self.metadata_reader.get_band_metadata_by_file_name(band_file)

        if not self.band_metadata:
            raise KeyError('Invalid band')


        self.band_dataset = gdal.Open(band_file)
        self.band_array = self.band_dataset.GetRasterBand(1).ReadAsArray()


    def get_radiance_as_array(self):
        radiance = ((self.band_metadata['radiance_maximum']-self.band_metadata['radiance_minimum']) / (self.band_metadata['quantize_cal_maximum']-self.band_metadata['quantize_cal_minimum'])) * (self.band_array - self.band_metadata['quantize_cal_minimum']) + self.band_metadata['radiance_minimum']
        radiance[self.band_array==0] = np.nan
        return radiance

    def get_reflectance_as_array(self, not_native_radiance_array=False):
        if self.band_metadata['type'] != 'reflectance':
            raise TypeError('Given band is thermal')
        if type(not_native_radiance_array)==bool:
            radiance = self.get_radiance_as_array()
        else:
            radiance = not_native_radiance_array
        d = float(self.metadata['EARTH_SUN_DISTANCE'])
        O = np.deg2rad(float(self.metadata['SUN_ELEVATION']))
        E = self.band_metadata['solar_irradiance']

        reflectance = (np.pi*radiance*d*d)/(E*np.sin(O))
        return reflectance

    def get_brightness_temperature_as_array(self):
        if self.band_metadata['type'] != 'thermal':
            raise TypeError('Given band is reflectance')

        radiance = self.get_radiance_as_array()

        K1 = self.band_metadata['k1_constant']
        K2 = self.band_metadata['k2_constant']

        brightness_temperature = (K2 / (np.log((K1/radiance+1))))

        return brightness_temperature


    def save_array_as_gtiff(self, array, new_file_path):
        driver = gdal.GetDriverByName("GTiff")
        dataType = gdal.GDT_Float32
        dataset = driver.Create(new_file_path, self.band_dataset.RasterXSize, self.band_dataset.RasterYSize, self.band_dataset.RasterCount, dataType)
        dataset.SetProjection(self.band_dataset.GetProjection())
        dataset.SetGeoTransform(self.band_dataset.GetGeoTransform())
        dataset.GetRasterBand(1).WriteArray(array)
        del dataset

class LandsatMetadataReader():

    def __init__(self, metadata_file_path):
        self.metadata_file_path = metadata_file_path
        self.metadata_file = open(self.metadata_file_path,'r')
        self.metadata = {}
        self.bands = {}

        solar_irradiances = {'LANDSAT_4':
                                 {'1':1957,
                                  '2':1825,
                                  '3':1557,
                                  '4':1033,
                                  '5':214.9,
                                  '7':80.72},
                             'LANDSAT_5':
                                 {'1': 1957,
                                  '2': 1826,
                                  '3': 1554,
                                  '4': 1036,
                                  '5': 215.0,
                                  '7': 80.67},
                             'LANDSAT_7':
                                 {'1': 1970,
                                  '2': 1842,
                                  '3': 1547,
                                  '4': 1044,
                                  '5': 225.7,
                                  '7': 82.06,
                                  '8': 1369},
                             'LANDSAT_8':
                                 {'1': 1895.33,
                                  '2': 2004.57,
                                  '3': 1820.75,
                                  '4': 1549.49,
                                  '5': 951.76,
                                  '6': 366.97,
                                  '7': 247.55,
                                  '8': 85.46,
                                  '9': 1723.88}
                             }

        wavelengths = {'LANDSAT_4':
                                {'1': 0.485,
                                 '2': 0.56,
                                 '3': 0.66,
                                 '4': 0.83,
                                 '5': 1.65,
                                 '6': 11.45,
                                 '7': 2.215},
                             'LANDSAT_5':
                                 {'1': 0.485,
                                  '2': 0.56,
                                  '3': 0.66,
                                  '4': 0.83,
                                  '5': 1.65,
                                  '6': 11.45,
                                  '7': 2.215},
                             'LANDSAT_7':
                                 {'1': 0.4775,
                                  '2': 0.560,
                                  '3': 0.6615,
                                  '4': 0.835,
                                  '5': 1.648,
                                  '6_VCID_1': 11.45,
                                  '6_VCID_2': 11.45,
                                  '7': 11.335,
                                  '8': 0.7055},
                             'LANDSAT_8':
                                 {'1': 0.44,
                                  '2': 0.48,
                                  '3': 0.56,
                                  '4': 0.655,
                                  '5': 0.865,
                                  '6': 1.61,
                                  '7': 2.20,
                                  '8': 0.59,
                                  '9': 1.37,
                                  '10': 10.9,
                                  '11': 12.05}
                        }

        center_sensor_zenith = {'LANDSAT_4':0.001,
                                'LANDSAT_5':0.001,
                                'LANDSAT_7':0.001,
                                'LANDSAT_8':0.001,}

        for line in self.metadata_file.readlines():
            if (line.find('GROUP') >= 0) or (line.find('=') == -1):
                continue
            else:
                line_normalized = line.replace(' ','')
                items = line_normalized.split('=')
                self.metadata[items[0]] = items[1].replace('\n','').replace('\"','')

        if not 'SPACECRAFT_ID' in self.metadata:
            raise KeyError('Invalid metadata file')

        all_bands = []
        reflectance_bands = []
        thermal_bands = []

        if self.metadata['SPACECRAFT_ID'] == 'LANDSAT_8':
            all_bands = [1,2,3,4,5,6,7,8,9,10,11]
            reflectance_bands = [1,2,3,4,5,6,7,8,9]
            thermal_bands = [10,11]

        if self.metadata['SPACECRAFT_ID'] == 'LANDSAT_7':
            all_bands = [1,2,3,4,5,'6_VCID_1','6_VCID_2',7,8]
            reflectance_bands = [1, 2, 3, 4, 5, 7, 8]
            thermal_bands = ['6_VCID_1','6_VCID_2']

        if self.metadata['SPACECRAFT_ID'] == 'LANDSAT_5':
            all_bands = [1, 2, 3, 4, 5, 6, 7]
            reflectance_bands = [1, 2, 3, 4, 5, 7]
            thermal_bands = [6]

        if self.metadata['SPACECRAFT_ID'] == 'LANDSAT_4':
            all_bands = [1, 2, 3, 4, 5, 6, 7]
            reflectance_bands = [1, 2, 3, 4, 5, 7]
            thermal_bands = [6]

        if not all_bands:
            raise KeyError('Invalid metadata file')

        for band in all_bands:
            self.bands[str(band)] = {}

            self.bands[str(band)]['file_name'] = self.metadata['FILE_NAME_BAND_%s' % band]
            self.bands[str(band)]['number'] = band
            self.bands[str(band)]['radiance_maximum'] = float(self.metadata['RADIANCE_MAXIMUM_BAND_%s' % band])
            self.bands[str(band)]['radiance_minimum'] = float(self.metadata['RADIANCE_MINIMUM_BAND_%s' % band])
            self.bands[str(band)]['quantize_cal_maximum'] = float(self.metadata['QUANTIZE_CAL_MAX_BAND_%s' % band])
            self.bands[str(band)]['quantize_cal_minimum'] = float(self.metadata['QUANTIZE_CAL_MIN_BAND_%s' % band])
            self.bands[str(band)]['radiance_mult'] = float(self.metadata['RADIANCE_MULT_BAND_%s' % band])
            self.bands[str(band)]['radiance_add'] = float(self.metadata['RADIANCE_ADD_BAND_%s' % band])

            self.bands[str(band)]['wavelength'] = wavelengths[self.metadata['SPACECRAFT_ID']][str(band)]

            self.bands[str(band)]['center_sensor_zenith'] = center_sensor_zenith[self.metadata['SPACECRAFT_ID']]


            if band in reflectance_bands:
                self.bands[str(band)]['saturation'] = self.metadata['SATURATION_BAND_%s' % band]
                self.bands[str(band)]['reflectance_maximum'] = float(self.metadata['REFLECTANCE_MAXIMUM_BAND_%s' % band])
                self.bands[str(band)]['reflectance_minimum'] = float(self.metadata['REFLECTANCE_MINIMUM_BAND_%s' % band])
                self.bands[str(band)]['reflectance_mult'] = float(self.metadata['REFLECTANCE_MULT_BAND_%s' % band])
                self.bands[str(band)]['reflectance_add'] = float(self.metadata['REFLECTANCE_ADD_BAND_%s' % band])

                self.bands[str(band)]['solar_irradiance'] = solar_irradiances[self.metadata['SPACECRAFT_ID']][str(band)]

                self.bands[str(band)]['type'] = 'reflectance'

            if band in thermal_bands:
                self.bands[str(band)]['k1_constant'] = float(self.metadata['K1_CONSTANT_BAND_%s' % band])
                self.bands[str(band)]['k2_constant'] = float(self.metadata['K2_CONSTANT_BAND_%s' % band])

                self.bands[str(band)]['type'] = 'thermal'

    def get_band_metadata_by_file_name(self, file_name):
        for band in self.bands.keys():
            if os.path.basename(file_name) == self.bands[band]['file_name']:
                return self.bands[band]

class LSTRetriever():
    c0 = 9.087
    c1 = 0.653
    c2 = -9.674
    cloudbufferdistance = 30
    cloudprobthreshold = 75.0
    shadowbufferdistance = 0
    NDVIs = 0.2
    NDVIv = 0.5

    LSE_modes = ['auto-ndvi-raw', 'auto-ndvi-srem', 'from-glc', 'external']

    def __init__(self, metadata_file, LSE_mode='auto-ndvi-srem', LSE_file=None, angles_file=None, usgs_utils=None,
                 temp_dir=None, window_size=7, cygwin_bash_exe_path=None):
        if LSE_mode not in self.LSE_modes:
            raise ValueError('Unsupported LSE mode. Supported: %s' % self.LSE_modes)

        self.metadata_file = metadata_file
        self.metadata = LandsatMetadataReader(self.metadata_file)

        self.dataset_basepath = os.path.dirname(self.metadata_file)
        self.window_size = window_size

        self.band_10_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_10'])
        self.band_11_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_11'])

        self.angles_file = angles_file
        self.usgs_utils = usgs_utils
        self.temp_dir = temp_dir
        self.lse_mode = LSE_mode
        self.cygwin_bash_exe_path = cygwin_bash_exe_path
        if self.lse_mode in ['from-glc', 'external']:
            self.lse_file = LSE_file

    def get_b_arrays(self, cwv_array):
        print('Creating b-coefs arrays')

        condition_list = [np.logical_and(cwv_array >= 0, cwv_array < 2.5),
                          np.logical_and(cwv_array >= 2.5, cwv_array < 3.5),
                          np.logical_and(cwv_array >= 3.5, cwv_array < 4.5),
                          np.logical_and(cwv_array >= 4.5, cwv_array < 5.5),
                          np.logical_and(cwv_array >= 5.5, cwv_array < 6.3),
                          cwv_array >= 6.3]

        choice_list_b0 = [-2.78009, 11.00824, 9.62610, 0.61258, -0.34808, -0.41165]
        choice_list_b1 = [1.01408, 0.95995, 0.96202, 0.99124, 0.98123, 1.00522]
        choice_list_b2 = [0.15833, 0.17243, 0.13834, 0.10051, 0.05599, 0.14543]
        choice_list_b3 = [-0.34991, -0.28852, -0.17262, -0.09664, -0.03518, -0.27297]
        choice_list_b4 = [4.04487, 7.11492, 7.87883, 7.85758, 11.96444, 4.06655]
        choice_list_b5 = [3.55414, 0.42684, 5.17910, 6.86626, 9.06710, -6.92512]
        choice_list_b6 = [-8.88394, -6.62025, -13.26611, -15.00742, -14.74085, -18.27461]
        choice_list_b7 = [0.09152, -0.06381, -0.07603, -0.01185, -0.20471, 0.24468]

        b0 = np.select(condition_list, choice_list_b0)
        b1 = np.select(condition_list, choice_list_b1)
        b2 = np.select(condition_list, choice_list_b2)
        b3 = np.select(condition_list, choice_list_b3)
        b4 = np.select(condition_list, choice_list_b4)
        b5 = np.select(condition_list, choice_list_b5)
        b6 = np.select(condition_list, choice_list_b6)
        b7 = np.select(condition_list, choice_list_b7)

        self.__save_array_to_gtiff(b0, gdal.Open(self.band_10_path), os.path.join(self.temp_dir, 'latest_bo.tif'))
        self.__save_array_to_gtiff(b3, gdal.Open(self.band_10_path), os.path.join(self.temp_dir, 'latest_b3.tif'))

        return b0, b1, b2, b3, b4, b5, b6, b7

    def get_lse_from_ndvi(self, srem=False):
        self.band_4_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_4'])
        self.band_5_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_5'])

        if not srem:
            band4_calibrator = CalibrateLandsatBand(self.band_4_path, self.metadata_file)
            band5_calibrator = CalibrateLandsatBand(self.band_5_path, self.metadata_file)

            band4_reflectance = band4_calibrator.get_reflectance_as_array()
            band5_reflectance = band5_calibrator.get_reflectance_as_array()

        else:
            from SREMPyLandsat.SREMPyLandsat import SREMPyLandsat

            srem = SREMPyLandsat(mode='landsat-usgs-utils')
            data_b4 = {'band': self.band_4_path,
                       'metadata': self.metadata_file,
                       'angles_file': self.angles_file,
                       'usgs_util_path': self.usgs_utils,
                       'temp_dir': self.temp_dir,
                       'cygwin_bash_exe_path': self.cygwin_bash_exe_path}

            srem.set_data(data_b4)
            band4_reflectance = srem.get_srem_surface_reflectance_as_array()

            data_b5 = {'band': self.band_5_path,
                       'metadata': self.metadata_file,
                       'angles_file': self.angles_file,
                       'usgs_util_path': self.usgs_utils,
                       'temp_dir': self.temp_dir,
                       'cygwin_bash_exe_path': self.cygwin_bash_exe_path}

            srem.set_data(data_b5)
            band5_reflectance = srem.get_srem_surface_reflectance_as_array()

        band4_reflectance[band4_reflectance > 1] = 1
        band5_reflectance[band5_reflectance > 1] = 1

        self.__save_array_to_gtiff(band4_reflectance, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_b4_reflectance.tif'))
        self.__save_array_to_gtiff(band5_reflectance, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_b5_reflectance.tif'))

        ndvi = (band5_reflectance - band4_reflectance) / (band5_reflectance + band4_reflectance)
        self.__save_array_to_gtiff(ndvi, gdal.Open(self.band_10_path), os.path.join(self.temp_dir, 'latest_ndvi.tif'))

        # from Land Surface Temperature Retrieval from Landsat 5, 7, and 8 over Rural Areas: Assessment of Different Retrieval Algorithms and Emissivity Models and Toolbox Implementation
        Pv = (ndvi - self.NDVIs) / (self.NDVIv - self.NDVIs)
        dE_b10 = (1 - 0.9668) * (1 - Pv) * 0.55 * 0.9863
        dE_b11 = (1 - 0.9747) * (1 - Pv) * 0.55 * 0.9896

        condition_list = [ndvi < self.NDVIs,
                          np.logical_and(ndvi >= self.NDVIs, ndvi <= self.NDVIv),
                          ndvi > self.NDVIv]

        low_ndvi_values_b10 = 0.973 - 0.047 * band4_reflectance
        medium_ndvi_values_b10 = 0.9863 * Pv + 0.9668 * (1 - Pv) + dE_b10
        high_ndvi_values_b10 = 0.9863 + dE_b10
        choice_list_b10 = [low_ndvi_values_b10, medium_ndvi_values_b10, high_ndvi_values_b10]

        low_ndvi_values_b11 = 0.984 - 0.0026 * band4_reflectance
        medium_ndvi_values_b11 = 0.9896 * Pv + 0.9747 * (1 - Pv) + dE_b11
        high_ndvi_values_b11 = 0.9896 + dE_b11
        choice_list_b11 = [low_ndvi_values_b11, medium_ndvi_values_b11, high_ndvi_values_b11]

        lse_b10 = np.select(condition_list, choice_list_b10)
        lse_b11 = np.select(condition_list, choice_list_b11)

        self.__save_array_to_gtiff(lse_b10, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_lse_b10.tif'))
        self.__save_array_to_gtiff(lse_b11, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_lse_b11.tif'))

        return lse_b10, lse_b11

    def lse_from_glc(self, glc_file_path):
        # Convert GLC to memory in domain
        base_raster = gdal.Open(self.band_10_path)
        geotransform = base_raster.GetGeoTransform()
        projection = base_raster.GetProjection()
        xMin = geotransform[0]
        yMax = geotransform[3]
        xMax = xMin + geotransform[1] * base_raster.RasterXSize
        yMin = yMax + geotransform[5] * base_raster.RasterYSize
        xRes = geotransform[1]
        yRes = geotransform[5]

        land_cover = gdal.Warp('', glc_file_path, format='MEM', xRes=xRes, yRes=yRes, dstSRS=projection,
                               outputBounds=[xMin, yMin, xMax, yMax])

        land_cover_array = land_cover.GetRasterBand(1).ReadAsArray()

        condition_list = [land_cover_array == 1,
                          land_cover_array == 2,
                          land_cover_array == 3,
                          land_cover_array == 4,
                          land_cover_array == 5,
                          land_cover_array == 6,
                          land_cover_array == 7,
                          land_cover_array == 8,
                          land_cover_array == 9,
                          land_cover_array == 10]

        choice_list_b10 = [0.971, 0.995, 0.97, 0.969, 0.992, 0.992, 0.98, 0.973, 0.969, 0.992]
        choice_list_b11 = [0.968, 0.996, 0.971, 0.97, 0.998, 0.998, 0.984, 0.981, 0.978, 0.998]

        lse_b10 = np.select(condition_list, choice_list_b10)
        lse_b11 = np.select(condition_list, choice_list_b11)

        self.__save_array_to_gtiff(lse_b10, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_lse_b10.tif'))
        self.__save_array_to_gtiff(lse_b11, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_lse_b11.tif'))

        return lse_b10, lse_b11

    def lse_from_file(self, lse_file_path):
        base_raster = gdal.Open(self.band_10_path)
        geotransform = base_raster.GetGeoTransform()
        projection = base_raster.GetProjection()
        xMin = geotransform[0]
        yMax = geotransform[3]
        xMax = xMin + geotransform[1] * base_raster.RasterXSize
        yMin = yMax + geotransform[5] * base_raster.RasterYSize
        xRes = geotransform[1]
        yRes = geotransform[5]

        lse = gdal.Warp('', lse_file_path, format='MEM', xRes=xRes, yRes=yRes, dstSRS=projection,
                               outputBounds=[xMin, yMin, xMax, yMax])

        lse_b10 = lse.GetRasterBand(0).ReadAsArray()
        lse_b11 = lse.GetRasterBand(1).ReadAsArray()

        return lse_b10, lse_b11


    def get_fmask_cloud_array(self):

        print('Running FMASK')
        output_cloud_path = os.path.join(self.temp_dir, 'latest_cloud.tif')
        cmd = '%s -o %s --scenedir %s --cloudbufferdistance %s --cloudprobthreshold %s --shadowbufferdistance %s' % (
        FMASK_EXECUTABLE_PATH, output_cloud_path, self.dataset_basepath, self.cloudbufferdistance, self.cloudprobthreshold,
        self.shadowbufferdistance)
        print('Command: %s' % cmd)
        os.system(cmd)
        cloud_ds = gdal.Open(output_cloud_path)
        cloud_array = cloud_ds.GetRasterBand(1).ReadAsArray()
        return cloud_array

    def get_cwv_as_array(self, cloud_array, b10_bt, b11_bt):

        b10_bt_masked = b10_bt.copy()
        b11_bt_masked = b11_bt.copy()
        print('Applying cloud/water masks')

        b10_bt_masked[cloud_array == 2] = np.nan
        b10_bt_masked[cloud_array == 5] = np.nan

        b11_bt_masked[cloud_array == 2] = np.nan
        b11_bt_masked[cloud_array == 5] = np.nan

        print('Splitting to subimages')
        b10_bt_subimages = self.get_subimages_with_step_and_size(b10_bt_masked, self.window_size, self.window_size)
        b11_bt_subimages = self.get_subimages_with_step_and_size(b11_bt_masked, self.window_size, self.window_size)

        print('Calculating CWV')
        cwv = self.compute_cwv_for_subimages(b10_bt_subimages, b11_bt_subimages)
        cwv[cwv < 0] = 0.0
        cwv[np.isnan(cwv)] = np.nanmean(cwv)

        print('Setting up CWV raster')
        cwv_full_array = self.transform_array_back_to_original_size(cwv, self.window_size, gdal.Open(self.band_10_path))

        cwv_full_array[cloud_array == 5] = np.nanmean(cwv_full_array)
        cwv_full_array[cwv_full_array < 0] = 0.0

        self.__save_array_to_gtiff(cwv_full_array, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_cwv.tif'))

        return cwv_full_array

    def get_subimages_with_step_and_size(self, image_array, step, size):
        subimages = []
        for r in range(0, image_array.shape[0] - size + 1, step):
            row_subimages = [image_array[r:r + size, c:c + size] for c in
                             range(0, image_array.shape[1] - size + 1, step)]
            subimages.append(row_subimages)

        return np.array(subimages)

    def compute_cwv_for_subimages(self, b10_subimages, b11_subimages):
        cwv_list = []
        i = 0
        print('Total subimages: %s' % len(b10_subimages))
        for row_b10, row_b11 in zip(b10_subimages, b11_subimages):
            cwv_row = []
            i += 1
            if i % 10 == 0:
                print('subimage: %s' % i)

            for item_b10, item_b11 in zip(row_b10, row_b11):
                cwp_item = self.compute_column_water_vapor_for_window(item_b10.reshape(-1), item_b11.reshape(-1))
                cwv_row.append(cwp_item)

            cwv_list.append(np.array(cwv_row))

        cwv_list = np.array(cwv_list)
        return cwv_list

    def get_b10_b11_brightness_temp_arrays(self):
        b10_calibrator = CalibrateLandsatBand(self.band_10_path, self.metadata_file)
        b11_calibrator = CalibrateLandsatBand(self.band_11_path, self.metadata_file)

        print('Band 10')
        b10_bt = b10_calibrator.get_brightness_temperature_as_array()
        print('Band 11')
        b11_bt = b11_calibrator.get_brightness_temperature_as_array()

        self.__save_array_to_gtiff(b10_bt, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_b10_bt.tif'))
        self.__save_array_to_gtiff(b11_bt, gdal.Open(self.band_10_path),
                                   os.path.join(self.temp_dir, 'latest_b11_bt.tif'))

        return b10_bt, b11_bt

    def get_lst_array(self):
        print('Start calculating LST')
        print('Calculating brightness temperatures')
        b10_bt, b11_bt = self.get_b10_b11_brightness_temp_arrays()

        print('Get clouds and water')
        cloud_array = self.get_fmask_cloud_array()
        cwv = self.get_cwv_as_array(cloud_array, b10_bt, b11_bt)
        b10_bt[cloud_array == 2] = np.nan

        # LSE
        if self.lse_mode == 'auto-ndvi-raw':
            lse_b10, lse_b11 = self.get_lse_from_ndvi(srem=False)
        if self.lse_mode == 'auto-ndvi-srem':
            lse_b10, lse_b11 = self.get_lse_from_ndvi(srem=True)
        if self.lse_mode == 'external':
            lse_b10, lse_b11 = self.lse_from_file(self.lse_file)
        if self.lse_mode == 'from-glc':
            lse_b10, lse_b11 = self.lse_from_glc(self.lse_file)

        lse_diff = lse_b10 - lse_b11
        lse_mean = (lse_b10 + lse_b11) / 2.0

        b0, b1, b2, b3, b4, b5, b6, b7 = self.get_b_arrays(cwv)

        print('Actually LST Split-Window formula')
        lst = b0 + (b1 + b2 * ((1 - lse_mean) / lse_mean) + b3 * (lse_diff / (lse_mean ** 2))) * (
                    (b10_bt + b11_bt) / 2.0) + (
                          b4 + b5 * ((1 - lse_mean) / lse_mean) + b6 * (lse_diff / (lse_mean ** 2))) * (
                          (b10_bt - b11_bt) / 2.0) + b7 * ((b10_bt - b11_bt) ** 2)
        return lst

    def get_lst_as_gtiff(self, output_path):
        lst = self.get_lst_array()
        self.__save_array_to_gtiff(lst, gdal.Open(self.band_10_path), os.path.join(self.temp_dir, output_path))

    def transform_array_back_to_original_size(self, array, window_size, original_raster):
        original_geotransform = original_raster.GetGeoTransform()
        original_projection = original_raster.GetProjection()
        original_datatype = gdal.GDT_Float32

        xMin = original_geotransform[0]
        yMax = original_geotransform[3]
        xMax = xMin + original_geotransform[1] * original_raster.RasterXSize
        yMin = yMax + original_geotransform[5] * original_raster.RasterYSize

        new_geotransform = (original_geotransform[0],
                            original_geotransform[1] * window_size,
                            original_geotransform[2],
                            original_geotransform[3],
                            original_geotransform[4],
                            original_geotransform[5] * window_size)

        driver = gdal.GetDriverByName('MEM')
        new_ds = driver.Create('', array.shape[1], array.shape[0], 1, original_datatype)
        new_ds.GetRasterBand(1).WriteArray(array)
        new_ds.SetProjection(original_projection)
        new_ds.SetGeoTransform(new_geotransform)

        # to source grid
        new_ds_translated = gdal.Warp('', new_ds, format='MEM', xRes=original_geotransform[1],
                                      yRes=original_geotransform[5], outputBounds=[xMin, yMin, xMax, yMax],
                                      dstNodata=np.nan, resampleAlg='cubic')

        new_ds_translated_array = new_ds_translated.GetRasterBand(1).ReadAsArray()

        return new_ds_translated_array

    def compute_column_water_vapor_for_window(self, tik, tjk):
        # feed with N pixels

        ti_mean = np.nansum(tik) / np.count_nonzero(~np.isnan(tik))
        tj_mean = np.nansum(tjk) / np.count_nonzero(~np.isnan(tjk))

        # numerator: sum of all (Tik - Ti_mean) * (Tjk - Tj_mean)
        numerator_ji_terms = []
        for ti, tj in zip(tik, tjk):
            numerator_ji_terms.append((ti - ti_mean) * (tj - tj_mean))

        numerator_ji = np.nansum(np.array(numerator_ji_terms)) * 1.0

        # denominator:  sum of all (Tik - Tj_mean)^2
        denominator_ji_terms = []
        for ti in tik:
            term = (ti - ti_mean) ** 2
            denominator_ji_terms.append(term)
        denominator_ji = np.nansum(np.array(denominator_ji_terms)) * 1.0

        # ratio ji
        ratio_ji = numerator_ji / denominator_ji

        # column water vapor
        cwv = self.c0 + self.c1 * (ratio_ji) + self.c2 * ((ratio_ji) ** 2)
        if math.isnan(float(cwv)):
            return np.nan

        return cwv

    def __save_array_to_gtiff(self, array, domain_raster, gtiff_path):
        driver = gdal.GetDriverByName("GTiff")
        dataType = gdal.GDT_Float32
        dataset = driver.Create(gtiff_path, array.shape[1], array.shape[0], domain_raster.RasterCount, dataType)
        dataset.SetProjection(domain_raster.GetProjection())
        dataset.SetGeoTransform(domain_raster.GetGeoTransform())

        dataset.GetRasterBand(1).WriteArray(array)

        del dataset
