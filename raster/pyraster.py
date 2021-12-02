# -*- coding: utf-8 -*-

import numpy as np
from numba import jit
from scipy import ndimage as ndi

import rasterio

@jit(nopython=True)
def calculate_slope(raster_input, np_easting, np_northing, noDataValue):
    """
    Calculates slope map based on an input raster map.
    
    Parameters
    ----------   
    raster_input (np.ndarray of float)
        Input raster (2D) where a slope map is calculated from
    np_easting (np.array of float)
        1D array of Easting coordinates for all pixels in the input raster, organised from left to right
    np_northing (np.array of float)
        1D array of Northing coordinates for all pixels in the input raster, organised from top to bottom
    noDataValue (float)
        Pixel value in the input raster corresponding to no data
    
    Return
    ------
    np.ndarray of float
        
    Assumed 3x3 window:

                        -------------------------
                        |   a   |   b   |   c   |
                        -------------------------
                        |   d   |   e   |   f   |
                        -------------------------
                        |   g   |   h   |   i   |
                        -------------------------
    """
    nRows, nCols = raster_input.shape;
    
    # absolute value of mean distance in meters between subsequent pixels in X & Y direction
    dx = np.abs((np_easting[1:] - np_easting[:-1]).mean());
    dy = np.abs((np_northing[1:] - np_northing[:-1]).mean());
    
    dzdx = np.full((nRows, nCols), 0, dtype=np.float32);
    dzdy = np.full((nRows, nCols), 0, dtype=np.float32);
    
    iter_northing = 0;
    for northing in np_northing:
        # first and last row of raster: skip
        if iter_northing == 0 or iter_northing == nRows-1:
            iter_northing += 1;
            continue;
        iter_easting = 0;
        for easting in np_easting:
            # first and last column of raster: skip
            # pixel value of raster corresponding to no data: skip
            if iter_easting == 0 or iter_easting == nCols-1\
                    or raster_input[iter_northing, iter_easting] == noDataValue:
                iter_easting += 1;
                continue;
            # moving 3x3 sub-raster
            a = raster_input[iter_northing-1, iter_easting-1]; # top left pixel value
            b = raster_input[iter_northing-1, iter_easting];
            c = raster_input[iter_northing-1, iter_easting+1]; # top right pixel value
            d = raster_input[iter_northing, iter_easting-1];
            e = raster_input[iter_northing, iter_easting]; # center pixel value
            f = raster_input[iter_northing, iter_easting+1];
            g = raster_input[iter_northing+1, iter_easting-1]; # bottom left pixel value
            h = raster_input[iter_northing+1, iter_easting];
            i = raster_input[iter_northing+1, iter_easting+1]; # bottom right pixel value
            
            # in case there is 1 pixel in moving 3x3 sub-raster corresponding to no data: skip
            if np.any(np.array([a,b,c,d,e,f,g,h,i]) == noDataValue):
                iter_easting += 1;
                continue;
            
            # gradient in X-direction
            dzdx[iter_northing, iter_easting] = ((c + 2*f + i)-(a + 2*d + g)) / (8*dx);
            # gradient in Y-direction
            dzdy[iter_northing, iter_easting] = ((g + 2*h + i)-(a + 2*b + c)) / (8*dy);
            iter_easting += 1;
        iter_northing += 1;
    hpot = np.hypot(np.abs(dzdy), np.abs(dzdx));
    return hpot;

def export_geotiff(raster_input, transform_out, filePath_out,
                   dtype_out=np.float32, epsg_out=31370, noDataValue=-9999):
    """
    Store a raster as geoTIFF 
    """
    # Transform raster data in specified data type
    raster = raster_input.astype(dtype_out);
    nRows, nCols = raster.shape;
    # Get a coordinate reference system based on the specified epsg number
    crs_out = rasterio.crs.CRS.from_epsg(epsg_out);
    bands=1;
    with rasterio.open(
            filePath_out, 'w', driver='GTiff',
            height=nRows, width=nCols,
            count=bands, dtype=dtype_out,
            crs=crs_out, transform=transform_out,
            nodata=noDataValue
    ) as outdata:
        for b in np.arange(bands):
            b += 1;
            b = np.int(b);
            if len(raster.shape) == 3:
                outdata.write(raster[:,:,b-1],b);
            else:
                outdata.write(raster,b);
        outdata.close();


class Topography(object):
    """Base class for general operations on 2D topography data"""
    def __init__(self, filePath_rasterTopo):
        """
        Parameters
        ----------
        filePath_rasterTopo (str)
            Full file path (incl. file extension) to location of 2D topography file.
        """
        self.filePath_rasterTopo = filePath_rasterTopo;
    
    def apply_rasterSlope(self):
        """Calculates 2D slope map from 2D topography map (raw data; not smoothed)."""
        rEasting_bottomRight = self.rEasting_topLeft\
                               + self.rPixelResolution_x*self.nCols;
        rNorthing_bottomRight = self.rNorthing_topLeft\
                                + self.rPixelResolution_y*self.nRows;
        np_easting = np.arange(
            self.rEasting_topLeft, rEasting_bottomRight, self.rPixelResolution_x
        );
        np_northing = np.arange(
            self.rNorthing_topLeft, rNorthing_bottomRight, self.rPixelResolution_y
        );
        self.rSlope = calculate_slope(
            self.rTopo, np_easting, np_northing, self.rNoDataValue
        );
    
    def get_rasterSlope(self, unit_out='Percentage'):
        """
        Returns the slope of the raster topography map (raw data; not smoothed). 
        
        Parameters
        ----------
        unit_out (str), optional
            Unit of the slope map {'Percentage', 'Radians', 'Degrees'}
            Default is 'Percentage'
        
        Returns
        -------
        numpy.ndarray (float)
        """
        if unit_out == 'Percentage':
            return 100*self.rSlope;
        elif unit_out == 'Radians':
            return np.arctan(self.rSlope);
        elif unit_out == 'Degree':
            return np.degrees(np.arctan(self.rSlope));

    def apply_rasterSlope2nd(self):
        """Calculates 2D second derivative map from 2D topography map (raw data; not smoothed)."""
        rEasting_bottomRight = self.rEasting_topLeft\
                               + self.rPixelResolution_x*self.nCols;
        rNorthing_bottomRight = self.rNorthing_topLeft\
                                + self.rPixelResolution_y*self.nRows;
        np_easting = np.arange(
            self.rEasting_topLeft, rEasting_bottomRight, self.rPixelResolution_x
        );
        np_northing = np.arange(
            self.rNorthing_topLeft, rNorthing_bottomRight, self.rPixelResolution_y
        );
        try:
            self.rSlope2nd = calculate_slope(
                self.rSlope, np_easting, np_northing, self.rNoDataValue
            );
        except AttributeError:
            # In case the slope was not yet calculated from the raster.
            self.apply_rasterSlope();
            self.rSlope2nd = calculate_slope(
                self.rSlope, np_easting, np_northing, self.rNoDataValue
            );
    
    def get_rasterSlope2nd(self):
        """
        Returns the second derivative map from 2D topography map (raw data; not smoothed)
            Unit = [100*m/m²]
        """
        return 100*self.rSlope2nd;

    def _rasterSmooth(self, raster_input, sigma_input=5.0, sigma_in_pxl=True):
        """See method apply_rasterSmooth"""
        if sigma_in_pxl:
            sigma = np.int(sigma_input);
        # If sigma_input is expressed in meter, a transformation to pixel is done first
        else:
            sigma_x = sigma_input / self.rPixelResolution_x;
            sigma_y = sigma_input / self.rPixelResolution_y;
            sigma = np.absolute([np.int(sigma_x), np.int(sigma_y)]);
        raster_input_smooth = ndi.gaussian_filter(
            raster_input, sigma, order=0, mode='reflect', truncate=2
        );
        return raster_input_smooth;
        
    def apply_rasterSmooth(
            self,
            rasterType='Topography',
            sigma_input=5.0,
            sigma_in_pxl=True,
            threshold_topo_low=0.5,
            threshold_topo_high=8.0):
        """
        Application of Gaussian smoothing on the selected raster.
        
        Parameters
        ----------
        rasterType (str), optional
            Type of raster: {'Topography', 'Slope', 'Slope2nd'}
            Default value is 'Topography'
        sigma_input (float or sequence of floats), optional
            Standard deviation for Gaussian kernel (single number = equal for all axes)
            Default value is 5.0
        sigma_in_pxl (bool), optional
            sigma_input value is expressed in pixels (True) or in meters (False)
            Default value is True
        threshold_topo_low (float or None), optional
            Replace raster values below threshold with no data value
             Only applicable when rasterType = 'Topography'
            Default value is 0.5
        threshold_topo_high (float or None), optional
            Replace raster values above threshold with no data value
             Only applicable when rasterType = 'Topography'
            Default value is 8.0
        """
        if rasterType == 'Topography':
            self.rTopo_smooth = self._rasterSmooth(
                self.rTopo, sigma_input, sigma_in_pxl
            );
            # Due to the high contrast in value between the topography values
            #  of the beach and the no data values (-9999),
            #  border effects of the smoothing take place at the transition.
            if threshold_topo_low is not None:
                self.rTopo_smooth[
                    self.rTopo_smooth <= threshold_topo_low] = self.rNoDataValue;
            if threshold_topo_high is not None:
                self.rTopo_smooth[
                    self.rTopo_smooth >= threshold_topo_high] = self.rNoDataValue;
        elif rasterType == 'Slope':
            self.rSlope_smooth = self._rasterSmooth(
                self.rSlope, sigma_input, sigma_in_pxl
            );
        elif rasterType == 'Slope2nd':
            self.rSlope2nd_smooth = self._rasterSmooth(
                self.rSlope2nd, sigma_input, sigma_in_pxl
            );
    
    def get_rasterRaw(self, rasterType='Topography', isMasked=False):
        """
        Return non-smoothed raster
        
        Parameters
        ----------
        rasterType (str), optional
            Type of raster: Topography (= DEFAULT), Slope
        isMasked (bool), optional
            output will be a masked array (False = DEFAULT)
                
        Returns
        -------
        numpy.ndarray (isMasked==FALSE) or numpy.ma.ndarray (isMasked==TRUE)
        """
        if rasterType == 'Topography' and isMasked:
            return self.rTopo;
        elif rasterType == 'Topography':
            return np.ma.masked_where(
                self.rTopo == self.rNoDataValue, self.rTopo, copy=True
            );
        elif rasterType == 'Slope' and isMasked:
            return self.rSlope;
        elif rasterType == 'Slope':
            return np.ma.masked_where(
                self.rSlope == self.rNoDataValue, self.rSlope, copy=True
            );

    def get_rasterSmooth(self, rasterType='Topography', isMasked=False):
        """
        Return smoothed raster
        
        Parameters
        ----------
        rasterType (str), optional
            Type of raster: Topography (= DEFAULT), Slope
        isMasked (bool), optional
            output will be a masked array (False = DEFAULT)
                
        Returns
        -------
        np.ndarray (isMasked==FALSE) or np.ma.ndarray (isMasked==TRUE)
        """
        if rasterType == 'Topography' and isMasked:
            return np.ma.masked_where(
                self.rTopo_smooth == self.rNoDataValue,
                self.rTopo_smooth, copy=True
            );
        elif rasterType == 'Topography':
            return self.rTopo_smooth;
        elif rasterType == 'Slope' and isMasked:
            return np.ma.masked_where(
                self.rSlope_smooth == self.rNoDataValue,
                self.rSlope_smooth, copy=True
            );
        elif rasterType == 'Slope':
            return self.rSlope_smooth;

    def apply_hinterlandFilter(self, threshold_topo=8.0):
        """
        Replace raster values above threshold with no data value
        
        Parameters
        ----------
        threshold_topo (float or None), optional
            Maximum threshold where smaller raster values are maintained
            Default value is 8.0
        """
        # If threshold_topo is None, no filtering is done in this method
        if threshold_topo is not None:
            self.threshTopo_hinterlandFilter = threshold_topo;
            self.rTopo[
                self.rTopo >= self.threshTopo_hinterlandFilter] = self.rNoDataValue;
            
    def apply_waterFilter(self, threshold_topo=0.5):
        """
        Replace raster values below threshold with no data value
        
        Parameters
        ----------
        threshTopo (float or None), optional
            Minimum threshold where higher raster values are maintained
            Default value is 0.5
        """
        # If threshold_topo is None, no filtering is done in this method
        if threshold_topo is not None:
            self.threshTopo_waterFilter = threshold_topo;
            self.rTopo[
                self.rTopo <= self.threshTopo_waterFilter] = self.rNoDataValue;

