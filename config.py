# -*- coding: utf-8 -*-

# Root directory where all data can be found : str
ROOT_DIRECTORY = r'';

# Type of platform in use : str
#  Examples: 'airplane', 'uav'
PLATFORM = 'airplane';
# Sensor used to make the acquisition : str
#  Examples: 'LiDAR', 'UX5HP_RGB15mm', 'EBEEX_MSREM'
SENSOR = 'LiDAR';
# Name of the site where the analysis should be done : str
#  Examples: 'Groenendijk', 'DeHaan', 'Mablethorpe', 'Dunkirk'
SITE_NAME = 'Groenendijk';


###############################################################################
RASTER_PREPROCESSING = {
    # File name of the shapefile : str or None
    #  File name of the shapefile used to crop the input raster (incl. file extension)
    #  If None, no cropping is performed on the input raster
    'CROP_ACCORDING_SHAPE_FILE': 'Groenendijk_LiDAR.shp',
    # Filter settings based on topography : dict of None or dict of float
    #  Lowest / highest topography value to take into account. Values lower/higher
    #   than this threshold are masked
    #  If None, no threshold is set
    'FILTER_BEACH_AREA_ON_RASTER': {
        'TOPOGRAPHY_THRESHOLD_LOW': 0.5, # expressed in meter
        'TOPOGRAPHY_THRESHOLD_HIGH': 8.0 # expressed in meter
    },
    'TOPOGRAPHY_GAUSSIAN_SMOOTHING': {
        # Should smoothing be done? : bool
        'DO': True,
        # Standard deviation of the Gaussian kernel : float
        'SIGMA': 2,
        # If the SIGMA parameter is defined in pixels (True) or meter (False) : bool
        'SIGMA_IN_PIXEL': False
    }
};


###############################################################################
CROSS_SECTION_SEARCH_BAR_FEATURES = {
    # Activation of the module? : bool
    'DO': True,
    # If the sea is located on the eastern side of the beach grid (e.g. east coast of UK),
    #  this parameter should be set to False, otherwise put True.
    'IS_SEA_AT_WESTERN_SIDE_OF_BEACH': True,
    # Orientation of the cross-sections relative to the beach, expressed in degree : float
    #  0    = Parallel to the beach
    #  90   = Perpendicular to the beach 
    'RELATIVE_ORIENTATION_DEGREE': 90,
    # Distance in meter between subsequent cross-section profile lines : float
    'INTER-SPACING_METER': 1,
    # Minimum length of valid data on the cross-section to take it into account: float
    'MINIMUM_LENGTH_CROSSSECTION_METER': 30, 
    # Minimum distance between crest and trough points in meter : float
    'MINIMUM_DISTANCE_CREST_TROUGH_METER': 3,
    # Minimum topography with lower values not to take into account : float or None
    #  If None, RASTER_PREPROCESSING['FILTER_BEACH_AREA_ON_RASTER']['TOPOGRAPHY_THRESHOLD_LOW'] is used
    #   If RASTER_PREPROCESSING['FILTER_BEACH_AREA_ON_RASTER']['TOPOGRAPHY_THRESHOLD_LOW'] is None, 'LOW_REFERENCE_LINE'[0] is used
    #    If 'LOW_REFERENCE_LINE'[0] is None, a default value of -5 is used
    'MINIMUM_TOPOGRAPHY_NOT_TO_TAKE_INTO_ACCOUNT_METER': None,
    # Frequency for plotting crest/trough points found on the cross-section profiles
    #  (expressed in number of cross-section profiles) : int
    # If -1, no plotting will be done
    'FREQUENCY_PLOTTING_CROSSSECTIONS': 50,
    # Topography intervals where to search for marker reference points : None or tuple of float
    #  If None, no searching is performed
    'LOW_REFERENCE_LINE': (0.6, 1.0),# expressed in meter
    'MEAN_LOW_WATER': None, # expressed in meter
    'MEAN_INTERTIDAL': None, # expressed in meter
    'MEAN_HIGH_WATER': None, # expressed in meter
    'HIGHEST_ASTRONOMICAL_TIDE': None # expressed in meter

};


###############################################################################
# Point filtering to remove isolated feature points on the raster
FILTERING_POINTS = {
    # Activation of the module? : bool
    'DO': False,
    # Topography feature point type where filter is applied on : str
    #  options: 'CrestPoints', 'TroughPoints' or 'InflectionPoints'
    'BEACH_FEATURE_TYPE': 'CrestPoints',
    # Dimensions of the filter window (width, height) in meter : tuple of float
    'FILTER_DIMENSIONS_METER': (20.0, 10.0),
    # Rotation angle (in degree) of the filtering window: None or float
    #  (+ in clockwise direction)
    #  if None, the main orientation of the beach is used (automatically determined)
    'FILTER_ROTATION_DEGREE': None,
    # Number of iterations the filter should be applied : int
    'NUMBER_ITERATIONS_TO_APPLY_FILTER_ON_RASTER': 3,
    # Number of points within filter window >= (k/100) * largest dimension of filter window
    'THRESHOLD_POINT_COUNT_IN_NEIGHBOURHOOD_PCT': 10, # =k [%] : int
    # Number of points within filter window with a topography difference > p (in meter)
    #  <= (q/100) * largest dimension of filter window
    'THRESHOLD_TOPOGRAPHY_DIFFERENCE_IN_NEIGHBOURHOOD': 0.10, # =p [m] : float
    'THRESHOLD_TOPOGRAPHY_DIFFERENCE_IN_NEIGHBOURHOOD_POINT_COUNT_PCT': 5 # =q [%] : int
};

        
###############################################################################
# Allocating a label number to each bar feature point on the raster, with points of the same
#  label thought to belong to the same individual intertidal bar / trough
GROUPING_POINTS_IN_LABELS = {
    # Activation of the module? : bool
    'DO': False,
    # Number of times the binary dilation should be repeated for respectively
    #  crest, trough and inflection points: list of int
    'NUMBER_ITERATIONS_DILATION_TO_GROUP_POINTS_METER': [10, 10, 3]
};
# Allocating a cluster label number to each bar feature point on the raster, with points of the
#  same cluster label thought to belong to the same cluster of intertidal bars / troughs
# The clustering is based on the DBSCAN algorithm (https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html)
#  and uses two features to cluster on:
#  1. Topography
#  2. Distance to a reference line parallel to the main direction of the intertidal bars
GROUPING_POINTS_IN_CLUSTERS = {
    # Activation of the module? : bool
    'DO': True,
    # Slope of the reference line to calculate the distance from : float or None
    #  Orientation (+) in clockwise direction
    #  If None, the main orientation of the beach is used
    'REFERENCE_LINE_SLOPE_DEGREE': None,
    # Point on the reference line to calculate the distance from : tuple of float
    #  Coordinates of the point (Easting, Northing) where the reference line
    #  should go through. Coordinate Reference System is the same as the raster
    'REFERENCE_LINE_POINT_COORDS': (55500, 219100),
    # Two main parameters of the DBSCAN clustering algorithm:
    #  First: eps
    #  Second: min_samples
    # For more information, please take a look at: scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    'DBSCAN_PARAMETERS': (0.10, 10),
    # Distance in topography to merge clusters, expressed in meter : float
    #  If two clusters are separated in topography (median value) within this distance,
    #  they are grouped together into one single cluster
    'MERGE_CLUSTERS_TOPOGRAPHY_INTERVAL': 0.3,
    # Percentage of total length of beach a cluster should be : float
    #  e.g. if parameter is 40 and a beach is 2000 meter in length,
    #  each cluster of points should be at least 800 meter in length
    'MINIMUM_SPATIAL_EXTENT_OF_CLUSTER_PCT': 0
}

                                    
###############################################################################
CROSS_SECTION_SEARCH_CHANNEL_FEATURES = {
    # bool
    'DO': True,
    # Orientation of the cross-sections relative to the beach, expressed in degree
    #  0    = Parallel to the beach
    #  90   = Perpendicular to the beach 
    'RELATIVE_ORIENTATION_DEGREE': 0,
    # Distance in meter between subsequent cross-section profile lines
    'INTER-SPACING_METER': 2,
    # Minimum length of valid data on the cross-section to take it into account: float
    'MINIMUM_LENGTH_CROSSSECTION_METER': 100,
    # Minimum width to consider a depression in topography as a channel, expressed in meter
    #  (measured between the top parts on both sides of the depression)
    'MINIMUM_WIDTH_METER': 3, #5 
    # Minimum depth to consider a depression in topography as a channel , expressed in meter
    #  (measured from the maximum to the minimum topography of the depression)
    'MINIMUM_DEPTH_METER': 0.2,
    'FREQUENCY_PLOTTING_CROSSSECTIONS': 10,
};