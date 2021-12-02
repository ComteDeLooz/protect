# -*- coding: utf-8 -*-

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
import rasterio as rio
import rasterio.mask
import rasterio.fill
import rasterio.crs
import fiona

import raster.pyraster
import constants

class DataAcquiredFromPlatform(raster.pyraster.Topography):
    """
    Specific methods for operations on LiDAR or UAV topography data.
    Inherits from the general Topography class.
    
    REQUIREMENT OF THE DATASET
        The length of the beach area (parallel to the water line) must be
         larger than the width of the beach area (perpendicular to the
         water line), otherwise the calculated beach orientation is not correct.
    """
    def __init__(self, filePath_topo, filePath_shapeFile=None):
        """
        Parameters
        ----------
        filePath_topo (str)
            Full file path (incl. file extension = .tif) to location of
             2D topography file from UAV or LiDAR measurements.
            No default is set.    
        filePath_shapeFile (str), optional
            Full file path (incl. file extension = .shp)) to location of
             shapefile defining the Region of Interest (ROI). 
             Default is None (no shapefile).
        """
        super().__init__(filePath_topo);
        self.read_topographyData(
            filePath_topo, filePath_shapeFile
        );
    
    def read_topographyData(self, filePath_raster, filePath_shapeFile=None):
        """
        Read LiDAR topography raster file (format: single-band geoTIFF), 
         potentially over a given area defined in a *.shp file.
        """
        with rio.open(filePath_raster) as dataset:
            # Both filePath_shapeFile and site_name need to be defined
            #   and valid to get the raster data over the ROI.
            if filePath_shapeFile is None:
                # Read the complete raster data
                self.rTopo = dataset.read(1, masked=False);
                self.rTransform_window = dataset.transform;
            else:
                # Read the raster data over a window specified by the shapefile
                with fiona.open(filePath_shapeFile, 'r') as shapefile:
                    shapes = [
                        feature['geometry'] for feature in shapefile
                    ];
                raster, self.rTransform_window = rasterio.mask.mask(
                    dataset, shapes, filled=True, crop=True
                );
                self.rTopo = raster[0,:,:];
            
            self.epsg_dataset = rasterio.crs.CRS(dataset.crs).to_epsg(); 
            self.rTransform_dataset = dataset.transform;
            self.rNoDataValue = dataset.nodatavals[0];
            
            self.rShape = self.rTopo.shape;
            self._rasterDimensions();
            self._pixelResolution();
            self._coordsTopLeft();
    
    def _rasterDimensions(self):
        self.nRows, self.nCols = self.rShape;   

    def _pixelResolution(self):
        """
        Extracts the pixel resolutions (NO absolute values) out of the
            affine transformation matrix for the window applied to raster.
        """
        self.rPixelResolution_x = self.rTransform_window[0];
        # Will be often negative due to the standard convention in image processing
        #   (top left corner of image is (0,0))
        self.rPixelResolution_y = self.rTransform_window[4]; 
    
    def _coordsTopLeft(self):
        """
        Extracts the coordinates of the top left corner of the dataset
            or shapefile window (if defined).
        """
        # If the shapefile window is completely within the dataset window,
        #   the top left coordinate of the shapefile window is taken.
        # If not: look for the top left coordinate the closest to the dataset.
        self.rEasting_topLeft = np.amax([
            self.rTransform_window[2], self.rTransform_dataset[2]
        ]);
        self.rNorthing_topLeft = np.amin([
            self.rTransform_window[5], self.rTransform_dataset[5]
        ]);
    
    def get_epsg(self):
        return self.epsg_dataset;
    
    def get_pixelResolution(self):
        """
        Returns
        -------
        float
            Number of meters corresponding to 1 pixel in X-direction
        float
            Number of meters corresponding to 1 pixel in Y-direction
            Is negative if Northing of top left pixel of raster topography
            is greater than Northing of bottom left pixel (almost always the case)
        """
        return self.rPixelResolution_x, self.rPixelResolution_y;
    
    def get_coordsRasterOrigin_refEN(self):
        """
        Returns
        -------
        float        
            Easting (CRS specified by geoTIFF) of maximum between top left pixel
            of raster topography and of shapefile window (if defined)
        float
            Northing (CRS specified by geoTIFF) of minimum between top left pixel
            of raster topography and of shapefile window (if defined)
        """
        return self.rEasting_topLeft, self.rNorthing_topLeft;
    
    def get_noDataValue(self):
        """
        Returns
        -------
        float
            No data value for the raster topography data
        """
        return self.rNoDataValue;
    
    def get_shape(self):
        """
        Returns
        -------
        tuple of int
            Dimensions of raster topography data (rows, columns)
            If valid window applies: raster dimensions of the cropped region
        """
        return self.rShape;
    
    def get_transformWindow(self):
        """
        Returns
        -------
        rio.Affine
            Affine transformation matrix for the window applied to raster
                based on shapefile specified as parameter
            If no window applies, this attribute will be the same as self.rTransform_dataset.
        """
        return self.rTransform_window;
    
    def apply_beachFilter(self):
        """
        This filtering follows on the earlier applied
         hinterland (filter out above threshold)
         and water filtering (filter out below threshold).
        The purpose is to maintain all pixels belonging to the beach
         and replace as much as possible all pixels not belonging
         to the beach with NoDataValue.
        """
        #######################################################################
        # Find the beach area in the raster data
        
        # Generate a binary raster (NoDataValue & 1)
        rTopo_binary = np.full(self.rShape, 1);
        rTopo_binary = np.where(
            self.rTopo != self.rNoDataValue,
            rTopo_binary,
            [self.rNoDataValue]
        );
        # Label all structures in the binary image as regions
        rTopo_label = label(rTopo_binary);
        # Calculate region properties
        rprops = regionprops(rTopo_label, intensity_image=rTopo_binary);
        area = 0;
        for rp in rprops:
            # Only store information when the largest region is found with values
            #  different from the no data value
            if rp.area > area and rp.mean_intensity != self.rNoDataValue:
                area = rp.area;
                label_beach = rp.label;
                bbox = rp.bbox;
                filled_image = rp.filled_image;
                image = rp.image;
        # Get all the pixels where the region label is equal to the largest found
        rTopo_beach = (rTopo_label == label_beach);
        # Set all pixels different from the largest region to no data value
        self.rTopo[~rTopo_beach] = self.rNoDataValue;
        # If there is a hole detected in the DEM
        if np.any(np.logical_xor(image, filled_image)):
            rTopo_sub = self.rTopo[
                bbox[0]:bbox[2],
                bbox[1]:bbox[3]
            ];
            rTopo_sub_mask = np.logical_and(
                (rTopo_sub == self.rNoDataValue),
                filled_image
            );
            self.rTopo[bbox[0]:bbox[2], bbox[1]:bbox[3]] = rasterio.fill.fillnodata(
                image=rTopo_sub,
                mask=np.logical_not(rTopo_sub_mask)
            );
        #######################################################################
        # Do binary closing on the selected beach area to remove "side effects"
        #  (e.g. artificial prolongation of the beach land inwards) and to keep
        #  only the main beach (better determination of beach orientation)
        
        # Re-generate a binary raster (NoDataValue & 1)
        rTopo_binary = np.full(self.rShape, 1);
        rTopo_binary = np.where(
            self.rTopo != self.rNoDataValue,
            rTopo_binary, [0]
        );
        # Binary closing over 100 iterations
        rTopo_binary_dilation = ndi.binary_dilation(
            rTopo_binary,
            iterations=constants.BEACH_FILTERING_BINARY_CLOSING_ITERATIONS
        );
        rTopo_binary_dilation_erosion = ndi.binary_erosion(
            rTopo_binary_dilation,
            iterations=constants.BEACH_FILTERING_BINARY_CLOSING_ITERATIONS
        );
        # Re-label all structures in the binary image as regions
        rTopo_label = label(rTopo_binary_dilation_erosion);
        # Calculate region properties
        rprops = regionprops(
            rTopo_label, intensity_image=rTopo_binary_dilation_erosion
        );
        area = 0;
        for rp in rprops:
            if rp.area > area and rp.mean_intensity != 0:
                # Based on the soft-bordered beach area, calculate orientation
                #  and length of the major axis
                orientation = rp.orientation
                self.major_axis_length = rp.major_axis_length;
        # Orientation in the function skimage.measure.regionprops is defined in the following axis system:
        #   X   0th axis (rows): applied to raster = vertically downward
        #   Y   1th axis (columns): applied to raster = horizontally to the right
        #   (+) in counter-clockwise direction
        #
        # In the PROTECT script, the axis system is defined differently:
        #   X   horizontally to the right
        #   Y   vertically downward
        #   (+) in clockwise direction
        #
        # Transforming angles from the first to the second axis system
        #   is obtained through: alpha' = (90° - alpha)
        self.beachOrientation = np.pi/2 - orientation;
        if self.beachOrientation < -np.pi/2:
            self.beachOrientation += np.pi;
        elif self.beachOrientation > np.pi/2:
            self.beachOrientation -= np.pi;

    def get_beachOrientation(self, unit='Radians'):
        """
        Angle of the main axis of the beach to the horizontal axis.
        Unit of angle is depending on argument value.
            +   in clockwise direction
            -   in counter-clockwise direction
        
        Arguments
        ---------
        unit (str), optional
            The unit of the orientation {'Radians', 'Degrees'}
        
        Returns
        -------
        float   
            Angle between the 0th axis (rows) and the major axis of the region,
                ranging from -pi/2 to pi/2 counter-clockwise.
        """
        if unit == 'Radians':
            return self.beachOrientation;
        elif unit == 'Degrees':
            return np.degrees(self.beachOrientation);

    def get_beachEstimatedLength(self):
        """Estimated length of the longest axis of the beach area"""
        return self.major_axis_length;
