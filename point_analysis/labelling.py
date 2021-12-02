# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import scipy.ndimage as ndi
import skimage.measure


class General(object):
    def __init__(self, filePath_pointCloud_csv, raster_shape):
        self.set_pointCloud(
            filePath_pointCloud_csv, raster_shape
        );

    def set_pointCloud(self, filePath_csv, raster_shape):
        # Read the feature point information from CSV-file
        self.pdPointCloud = pd.read_csv(
            filePath_csv, sep=';', index_col=0
        );
        # Generate a raster with values zero except for the locations of the
        #  feature points
        self.rPointCloud = np.full(
            raster_shape, 0, dtype=np.bool_
        );
        point_xy = self.pdPointCloud.loc[:,['Point_X','Point_Y']].values;
        point_xy = point_xy.astype(np.int);
        self.rPointCloud[point_xy[:,1], point_xy[:,0]] = 1;


class FromPointToGroupLabel(General):
    def __init__(
            self, filePath_pointCloud_csv, raster_shape, pixel_resolution_median):
        super().__init__(
            filePath_pointCloud_csv, raster_shape
        );
        self.pixel_resolution_median = pixel_resolution_median;
                
    def group_points_in_labels(self, number_iterations=3):
        # Apply a binary dilation to the points in the raster
        rPointCloud_dil = ndi.binary_dilation(
            self.rPointCloud, iterations=number_iterations
        );
        # Unique label for all groups of points
        labels = skimage.measure.label(
            rPointCloud_dil, return_num=False
        );
        # Paste on all original non-zero elements of the raster the label number
        self.rPointCloud_labels = np.multiply(self.rPointCloud, labels);
        # Add the labels for each non-zero element of the raster to the Pandas DF
        npLabel = np.array([]);
        pdPoint_xy = self.pdPointCloud.loc[:,['Point_X','Point_Y']];
        for pdIx, row in pdPoint_xy.iterrows():
            pointLabel = self.rPointCloud_labels[
                np.int(row.loc['Point_Y']), np.int(row.loc['Point_X'])
            ];
            npLabel = np.append(npLabel, pointLabel);
        self.pdPointCloud.loc[:,'Group_Label'] = npLabel;
    
    def calc_intertidalBarCharacteristics(self):
        self.pdPointCloud.loc[:,'Intertidal_Bar_Centroid_X'] = -1;
        self.pdPointCloud.loc[:,'Intertidal_Bar_Centroid_Y'] = -1;
        self.pdPointCloud.loc[:,'Intertidal_Bar_Orientation_degrees'] = -1;
        self.pdPointCloud.loc[:,'Intertidal_Bar_Width_m'] = -1;
        
        regionprops = skimage.measure.regionprops(self.rPointCloud_labels);
        for rp in regionprops:
            npBarLabel = rp.label;
            self.pdPointCloud.loc[
                self.pdPointCloud['Group_Label'] == npBarLabel,
                'Intertidal_Bar_Centroid_X'
            ] = np.int(rp.centroid[1]);
            self.pdPointCloud.loc[
                self.pdPointCloud['Group_Label'] == npBarLabel,
                'Intertidal_Bar_Centroid_Y'
            ] = np.int(rp.centroid[0]);
            # For difference in angles: see method apply_beachFilter of class DataAcquiredFromPlatform
            bar_orientation = np.pi/2 - rp.orientation;
            if bar_orientation < -np.pi/2:
                bar_orientation += np.pi;
            elif bar_orientation > np.pi/2:
                bar_orientation -= np.pi;
            self.pdPointCloud.loc[
                self.pdPointCloud['Group_Label'] == npBarLabel,
                'Intertidal_Bar_Orientation_degrees'
            ] = np.around(np.degrees(bar_orientation), 1);
            self.pdPointCloud.loc[
                self.pdPointCloud['Group_Label'] == npBarLabel,
                'Intertidal_Bar_Width_m'
            ] = np.int(rp.major_axis_length * self.pixel_resolution_median);
    
    def get_appendedInformation(self):
        return self.pdPointCloud;


class BottomPointToChannel(General):
    def __init__(self, filePath_pointCloud_csv, processingObject):
        super().__init__(filePath_pointCloud_csv, processingObject);





