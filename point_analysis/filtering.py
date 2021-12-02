# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from skimage.draw import polygon

def distance(point1, point2):
    return np.sqrt(
        np.power(point1[0]-point2[0], 2) + np.power(point1[1]-point2[1], 2)
    );


class NeighbourhoodFilter(object):
    """
    The main focus of the class is to filter out isolated points
    
    Parameters
    ----------
        filePath_pointCloud_csv (str)
            Full file path (incl. .csv) of the csv-file containing
                the location of identified points (crest or trough or inflection point)
                together with other information like topography, point ID or cross-section ID 
            The following columns are mandatory:
                * Point ID
                * Position in Easting, Northing
                * Position in X, Y
                * Topography at the point location
        raster_pixel_resolution (tuple of float)
            Pixel resolution in x-direction [0] and y-direction [1]
        raster_shape (tuple of int)
            Shape of the raster (rows, columns)
    """
    def __init__(
            self,
            filePath_pointCloud_csv,
            raster_pixel_resolution,
            raster_shape):
        
        rPixelResolution_x, rPixelResolution_y = raster_pixel_resolution;
        # Taking absolute values because most often rPixelResolution_y will be
        #   negative (origin image coordinates is top left pixel)
        self.rPixelResolution_mean = np.mean(
            np.absolute([rPixelResolution_x, rPixelResolution_y])
        );
        self.rShape = raster_shape;
        
        # Read csv-file and identify the columns of importance in the class.
        pdPointCloudInfo = pd.read_csv(
            filePath_pointCloud_csv, sep=';', index_col=0
        );
        self.colsHeader_original = pdPointCloudInfo.columns;
        # Find the indices of several column headers
        self.ixCol_pointID = np.argwhere(
            self.colsHeader_original == 'Point_ID'
        )[0][0];
        self.ixCol_easting = np.argwhere(
            self.colsHeader_original == 'Point_Easting'
        )[0][0];
        self.ixCol_northing = np.argwhere(
            self.colsHeader_original == 'Point_Northing'
        )[0][0];
        self.ixCol_x = np.argwhere(
            self.colsHeader_original == 'Point_X'
        )[0][0];
        self.ixCol_y = np.argwhere(
            self.colsHeader_original == 'Point_Y'
        )[0][0];
        self.ixCol_topo = np.argwhere(
            self.colsHeader_original == 'Point_Topography'
        )[0][0];
        # Transform Pandas Dataframe to Numpy array
        self.npPointCloudInfo = pdPointCloudInfo.to_numpy();
        # Creation of a zero-valued raster with the same shape as the topography raster
        self.rPointCloudTopo = np.zeros(self.rShape, dtype=np.float32);
        # All elements in the new raster remain zero except at the point locations
        #   listed in the csv-file given as parameter (feature points).
        # The values of the raster are the topography values at the feature point
        #   locations also listed in the same csv-file.
        self.rPointCloudTopo[
            self.npPointCloudInfo[:,self.ixCol_y].astype(np.int),
            self.npPointCloudInfo[:,self.ixCol_x].astype(np.int)
        ] = self.npPointCloudInfo[:,self.ixCol_topo];
        
        self._init_table_with_filter_info();
        
    def _init_table_with_filter_info(self):
        # List containing statistical information of the neighbourhood filtering
        #  per original point provided as input (thus also the points filtered out)
        self.npStats_neighbourhood = np.full(
            (self.npPointCloudInfo.shape[0], 19), -1, dtype=np.float32
        );
        # Copy information from input csv in new array
        self.npStats_neighbourhood[:,:4] = self.npPointCloudInfo[:, [
            self.ixCol_pointID,
            self.ixCol_easting,
            self.ixCol_northing, 
            self.ixCol_topo
        ]];
        
    def set_filter_characteristics(
            self, neighbourhoodDimensions_m=(20,10), rotationAngle_deg=0,
            number_iterations=3, threshold_point_count=0.3,
            threshold_topo_difference=0.10,
            threshold_topo_difference_point_count=0.15):
        """
        Add topography information
        Criteria to be fulfilled to withheld the target point (not filtered out)
            * neighbourhoodThreshold_countNonzero = 0.3 * self.neighbourhoodDimensions_pxl[0]
                The number of other points within the defined rectangular neighbourhood
                    should be at least 30% of the horizontal dimension expressed in pixels of the rectangular neighbourhood (x-axis)
            * neighbourhoodThreshold_diffTopo = 0.10
              neighbourhoodThreshold_countDiffTopo = 0.15 * self.neighbourhoodDimensions_pxl[0]
                The number of other points within the defined rectangular neighbourhood & with a difference in topography with the target point greater than neighbourhoodThreshold_diffTopo
                    should be less than 15% of the horizontal dimension expressed in pixels of the rectangular neighbourhood (x-axis)
                This second criterion is only tested during the first iteration (if more iterations apply)
        """
        self.neighhourhoodDimensions_m = neighbourhoodDimensions_m;
        # Transform neighbourhood dimensions from meter to pixel
        self.neighbourhoodDimensions_pxl = (
            neighbourhoodDimensions_m[0] / self.rPixelResolution_mean,
            neighbourhoodDimensions_m[1] / self.rPixelResolution_mean
        );
        # Transform rotation angle from degrees to radians
        self.rotationAngle_rad = np.radians(rotationAngle_deg);
        self.number_iterations = number_iterations;
    
        # Set-up of criteria to decide if the target point is filtered out
        #   Minimum number of points within the defined horizontal neighbourhood
        #    of the target point
        #   Defined relative to the size of the filter in x-direction
        self.neighbourhood_countNonzero_threshold = threshold_point_count\
                                                    * self.neighbourhoodDimensions_pxl[0];
        #   Maximum difference in topography between target point and neighbouring points
        #    This is the definition of neighbourhood in the Z-direction
        self.neighbourhood_diffTopo_threshold = threshold_topo_difference;
        #   Minimum number of points within the defined vertical neighbourhood
        #    of the target point
        #   Defined relative to the size of the filter in x-direction
        self.neighbourhood_countDiffTopo_threshold = threshold_topo_difference_point_count\
                                                     *self.neighbourhoodDimensions_pxl[0];
        # Store the filter characteristics in the statistical overview file
        self.npStats_neighbourhood[:,4:8] = [
            number_iterations,
            neighbourhoodDimensions_m[0], neighbourhoodDimensions_m[1],
            rotationAngle_deg
        ];
        # Add the criteria to the list containing statistical information of the
        #  neighbourhood filtering
        self.npStats_neighbourhood[:,[9,11,13]] = [
            self.neighbourhood_countNonzero_threshold,
            self.neighbourhood_countDiffTopo_threshold,
            self.neighbourhood_diffTopo_threshold
        ];
    
    def _rectangular_neighbourhood(self, pointCoords_refImage):
        """
        Sets pixel coordinates of the rotated rectangle with user-specified size
         (defined as a polygon).
        
        Arguments
        ---------
        pointCoords_refImage (tuple of int)
            Image coordinates in pixels (x,y) of the target point where a rotated
             rectanglular neighbourhood needs to be taken into account.
            
            X   horizontally to the right
            Y   vertically downward
            (+) in clockwise direction
        """
        size_x, size_y = self.neighbourhoodDimensions_pxl;
        # Coordinates of the non-rotated rectangle defined in the reference
        #  system with origin in the target point (expressed in pixels)
        rectangularCorners_refRec_noRot = [
            np.array([size_x/2, size_y/2]),
            np.array([size_x/2, -size_y/2]),
            np.array([-size_x/2, -size_y/2]),
            np.array([-size_x/2, size_y/2])
        ];
        # Define rotation matrix based on specified angles
        rotationMatrix = np.array([
            [np.cos(self.rotationAngle_rad), -np.sin(self.rotationAngle_rad)],
            [np.sin(self.rotationAngle_rad), np.cos(self.rotationAngle_rad)]
        ], dtype=np.float32);   
        
        rectangularCorners_refRec_rot = np.zeros(
            (len(rectangularCorners_refRec_noRot), 2)
        );
        # Apply rotation to the four corners of the rectangle. Coordinates are
        #  still defined in the reference system with origin in the target point
        for row_count, cornerCoords in enumerate(rectangularCorners_refRec_noRot):
            rectangularCorners_refRec_rot[row_count, :] = np.matmul(
                rotationMatrix, cornerCoords
            );
        
        rectangularCorners_refImage = np.zeros(
            (len(rectangularCorners_refRec_noRot), 2),
            dtype=np.int
        );
        # Translate coordinates in the reference system with origin in the target
        #  point to the reference system of the image (origin at top left). This
        #  is done by a translation operation.
        rectangularCorners_refImage[:,0] = pointCoords_refImage[0]\
                                           + rectangularCorners_refRec_rot[:,0];
        rectangularCorners_refImage[:,1] = pointCoords_refImage[1]\
                                           + rectangularCorners_refRec_rot[:,1];
        # Generate a polygon through skimage.draw
        self.neighbourhoodPolygon_rr, self.neighbourhoodPolygon_cc = polygon(
            rectangularCorners_refImage[:,1], rectangularCorners_refImage[:,0],
            shape=self.rShape
        );

    def apply_rectangular_filter(self):
        # Take a copy of the original numpy array because rows will be deleted
        #  and the original version is maintained. 
        self.npPointCloudInfo_filtered = np.copy(self.npPointCloudInfo);
        
        # Initialisation of a numpy array containing the indices to be deleted.
        #  This numpy array is re-initialised after each iteration in the loop.
        index_toBeDeleted = np.array([], dtype=np.uint);
        for it in np.arange(self.number_iterations):            
            count_it = 0;
            for npPoint in self.npPointCloudInfo_filtered:
                isPointFilteredOut = False;
                
                point_ID = np.int(npPoint[self.ixCol_pointID]);
                point_x = np.int(npPoint[self.ixCol_x]); 
                point_y = np.int(npPoint[self.ixCol_y]);
                
                self._rectangular_neighbourhood((point_x, point_y));
                # Get the values (1D array) of self.rPointCloudTopo in the
                #  rotated rectangular neighbourhood of the target point
                npPointCloudTopo_neighbourhood = self.rPointCloudTopo[
                    self.neighbourhoodPolygon_rr, self.neighbourhoodPolygon_cc
                ];
                # Count the number of pixels with a value different from zero
                neighbourhood_countNonzero = np.count_nonzero(
                    npPointCloudTopo_neighbourhood
                );
                # If the number of non-zero elements is less than a certain value
                #  (see method self.set_filter_characteristics):
                #  the point is considered to be isolated.
                if neighbourhood_countNonzero < self.neighbourhood_countNonzero_threshold:
                    # To filter the point out, its pixel value is set to zero
                    self.rPointCloudTopo[point_y, point_x] = 0;
                    index_toBeDeleted = np.append(index_toBeDeleted, count_it);
                    isPointFilteredOut = True;
                elif it == 0:
                    # Maintain only the useful topography values (!= 0) 
                    npPointCloudTopo_neighbourhood_nonzero = npPointCloudTopo_neighbourhood[
                        npPointCloudTopo_neighbourhood != 0
                    ];
                    # Topography at the location of the target point
                    topo_target = npPoint[self.ixCol_topo];
                    
                    # Identify all points in the rectangular neighbourhood around
                    #  the target point that have a topography difference with
                    #  the point more than a certain value (see method self.set_filter_characteristics)
                    npPointCloudTopo_neighbourhood_threshold = np.absolute(
                        npPointCloudTopo_neighbourhood_nonzero - topo_target
                    ) > self.neighbourhood_diffTopo_threshold;
                    pointCloudTopo_neighbourhood_thresholdCount = np.sum(
                        npPointCloudTopo_neighbourhood_threshold
                    );
                    # If the number of points in the neighbourhood with a too large
                    #  topography difference distance, is greater then a threshold,
                    #  the target point is filtered out
                    if pointCloudTopo_neighbourhood_thresholdCount > self.neighbourhood_countDiffTopo_threshold:
                        self.rPointCloudTopo[point_y, point_x] = 0;
                        index_toBeDeleted = np.append(index_toBeDeleted, count_it);
                        isPointFilteredOut = True;
                
                # Execute only in the last iteration or when the target point
                #  is filtered out
                if (it == self.number_iterations-1) or isPointFilteredOut:
                    npFindRow = self.npStats_neighbourhood[:,0] == point_ID;
                    self.npStats_neighbourhood[npFindRow, 8] = isPointFilteredOut;
                    
                    # Maintain only the useful topography values (!= 0) 
                    npPointCloudTopo_neighbourhood_nonzero = npPointCloudTopo_neighbourhood[
                        npPointCloudTopo_neighbourhood != 0
                    ];
                    # Topography at the location of the target point
                    topo_target = npPoint[self.ixCol_topo];
                    # Identify all points in the rectangular neighbourhood around
                    #  the target point that have a topography difference with
                    #  the point more than a certain value (see method self.set_filter_characteristics)
                    npPointCloudTopo_neighbourhood_threshold = np.absolute(
                        npPointCloudTopo_neighbourhood_nonzero - topo_target
                    ) > self.neighbourhood_diffTopo_threshold;
                    pointCloudTopo_neighbourhood_thresholdCount = np.sum(
                        npPointCloudTopo_neighbourhood_threshold
                    );
                    # Store certains countings over the neighbourhood filter
                    self.npStats_neighbourhood[npFindRow,[10,12]] = np.array([
                        neighbourhood_countNonzero,
                        pointCloudTopo_neighbourhood_thresholdCount
                    ]);
                    self.npStats_neighbourhood[npFindRow,14:] = np.around(np.array([
                        np.mean(pointCloudTopo_neighbourhood_thresholdCount),
                        np.median(pointCloudTopo_neighbourhood_thresholdCount),
                        np.amin(pointCloudTopo_neighbourhood_thresholdCount),
                        np.amax(pointCloudTopo_neighbourhood_thresholdCount),
                        np.std(pointCloudTopo_neighbourhood_thresholdCount)
                    ], dtype=np.float32), 4);
                count_it += 1;
            # After each iteration, delete the target points found to be filtered out
            self.npPointCloudInfo_filtered = np.delete(
                self.npPointCloudInfo_filtered, index_toBeDeleted, 0
            );
            # Re-initialize the array after each iteration
            index_toBeDeleted = np.array([], dtype=np.uint);

    def get_filteredPointsInTableForm(self):
        return pd.DataFrame(
            self.npPointCloudInfo_filtered, columns=self.colsHeader_original
        );
        
    def get_filterInfoTable(self):
        colsHeader = np.array([
            'Point_ID', 'Point_Easting', 'Point_Northing', 'Point_Topography',
            'Number_FilterIterations', 'NeighbourhoodDims_x_m', 'NeighbourhoodDims_y_m',
            'NeighbourhoodRotationAngle_deg', 'IsPointFilteredOut',
            'Threshold_Count_NeighbourhoodPoints', 'Count_NeighbourhoodPoints',
            'Threshold_Count_TopographyDifference', 'Count_TopographyDifference',
            'Threshold_Max_TopographyDifference', 'Mean_TopographyDifference',
            'Median_TopographyDifference', 'Min_TopographyDifference',
            'Max_TopographyDifference', 'STD_TopographyDifference'
        ]); 
        return pd.DataFrame(self.npStats_neighbourhood, columns=colsHeader);







