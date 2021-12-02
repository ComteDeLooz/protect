# -*- coding: utf-8 -*-

import os
from osgeo import gdal
import rasterio.crs

def subsample(
    filePath_raster_input, filePath_raster_output,
    target_resolution_x, target_resolution_y):
    if os.path.exists(filePath_raster_input):
        print('Input File exists');

    ds = gdal.Open(filePath_raster_input);
    gdal.Translate(
        filePath_raster_output,
        ds,
        format='GTiff',
        xRes=target_resolution_x,
        yRes=target_resolution_y,
        noData=-9999,
        resampleAlg='bilinear'
    );

def convert_adf_to_geotiff(
        filePath_adf, filePath_geotiff, epsg_number):
    ds = gdal.Open(filePath_adf);
    crs = rasterio.crs.CRS.from_epsg(epsg_number);
    gdal.Translate(
        filePath_geotiff,
        ds,
        format='GTiff',
        outputSRS=crs,
        noData=-9999
    );
