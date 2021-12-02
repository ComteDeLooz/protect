# -*- coding: utf-8 -*-

import logging
import os
import config

logger = logging.getLogger(__name__);
logger.info('Info log level - project.py');
logger.debug('Debug log level - project.py');
logger.critical('Critical log level - project.py');

class User(object):
    def __init__(self):
        self._init_raster_preprocessing();
        self._init_crossSection_barFeatures();
        self._init_crossSection_channelFeatures();
        self._init_pointFiltering();
        self._init_point_grouping_into_labels();
        self._init_point_clustering();
    
    ###########################################################################
    # USER INPUT RELATED TO RASTER PRE-PROCESSING
    ###########################################################################
    def _init_raster_preprocessing(self):
        preprocessing = config.RASTER_PREPROCESSING;

        self.set_filterBeachAreaOnRaster_topoThresholdLow(
            preprocessing['FILTER_BEACH_AREA_ON_RASTER']['TOPOGRAPHY_THRESHOLD_LOW']
        );
        self.set_filterBeachAreaOnRaster_topoThresholdHigh(
            preprocessing['FILTER_BEACH_AREA_ON_RASTER']['TOPOGRAPHY_THRESHOLD_HIGH']
        );
        self.set_isDoTopoSmoothing(
            preprocessing['TOPOGRAPHY_GAUSSIAN_SMOOTHING']['DO']
        );
        self.set_topoSmoothingSigma(
            preprocessing['TOPOGRAPHY_GAUSSIAN_SMOOTHING']['SIGMA']
        );
        self.set_isTopoSmoothingSigmaExpressedInPixels(
            preprocessing['TOPOGRAPHY_GAUSSIAN_SMOOTHING']['SIGMA_IN_PIXEL']
        );
    
    def set_filterBeachAreaOnRaster_topoThresholdLow(self, float_):
        self.filterBeachAreaOnRaster_topoThresholdLow = float_;
    
    def get_filterBeachAreaOnRaster_topoThresholdLow(self):
        return self.filterBeachAreaOnRaster_topoThresholdLow;
    
    def set_filterBeachAreaOnRaster_topoThresholdHigh(self, float_):
        self.filterBeachAreaOnRaster_topoThresholdHigh = float_;
    
    def get_filterBeachAreaOnRaster_topoThresholdHigh(self):
        return self.filterBeachAreaOnRaster_topoThresholdHigh;
    
    def set_isDoTopoSmoothing(self, flag):
        self.isDoTopoSmoothing = flag;
    def get_isDoTopoSmoothing(self):
        return self.isDoTopoSmoothing;
    
    def set_topoSmoothingSigma(self, float_):
        self.topoSmoothingSigma = float_;
    
    def get_topoSmoothingSigma(self):
        return self.topoSmoothingSigma;
    
    def set_isTopoSmoothingSigmaExpressedInPixels(self, flag):
        self.isTopoSmoothingSigmaExpressedInPixels = flag;
    
    def get_isTopoSmoothingSigmaExpressedInPixels(self):
        return self.isTopoSmoothingSigmaExpressedInPixels;
    
    ###########################################################################
    # USER INPUT RELATED TO INTERTIDAL BAR FEATURES
    ###########################################################################
    def _init_crossSection_barFeatures(self):
        crossSection_barFeatures = config.CROSS_SECTION_SEARCH_BAR_FEATURES
        self.set_isDoCrossSectionSearchBarFeatures(
            crossSection_barFeatures['DO']
        );
        self.set_crossSection_barFeatures_isSeaAtWesternSideOfBeach(
            crossSection_barFeatures['IS_SEA_AT_WESTERN_SIDE_OF_BEACH']
        );
        self.set_crossSection_barFeatures_relativeOrientation_deg(
            crossSection_barFeatures['RELATIVE_ORIENTATION_DEGREE']
        );
        self.set_crossSection_barFeatures_interSpacing_meter(
            crossSection_barFeatures['INTER-SPACING_METER']
        );
        self.set_crossSection_barFeatures_min_length_crossSection_meter(
            crossSection_barFeatures['MINIMUM_LENGTH_CROSSSECTION_METER']
        );
        self.set_crossSection_barFeatures_min_distance_crest_trough_meter(
            crossSection_barFeatures['MINIMUM_DISTANCE_CREST_TROUGH_METER']
        );
        self.set_crossSection_barFeatures_min_topo_not_take_into_account_meter(
            crossSection_barFeatures['MINIMUM_TOPOGRAPHY_NOT_TO_TAKE_INTO_ACCOUNT_METER']
        );
        self.set_searchIntervalLowReferenceLine(
            crossSection_barFeatures['LOW_REFERENCE_LINE']
        );
        self.set_searchIntervalMeanLowWater(
            crossSection_barFeatures['MEAN_LOW_WATER']
        );
        self.set_searchIntervalMeanIntertidal(
            crossSection_barFeatures['MEAN_INTERTIDAL']
        );
        self.set_searchIntervalMeanHighWater(
            crossSection_barFeatures['MEAN_HIGH_WATER']
        );
        self.set_searchIntervalHighestAstronomicalTide(
            crossSection_barFeatures['HIGHEST_ASTRONOMICAL_TIDE']
        );
        self.set_crossSection_plotting_frequency(
            crossSection_barFeatures['FREQUENCY_PLOTTING_CROSSSECTIONS']
        );
    
    def set_isDoCrossSectionSearchBarFeatures(self, flag):
        self.isDoCrossSectionSearchBarFeatures = flag;
    
    def get_isDoCrossSectionSearchBarFeatures(self):
        return self.isDoCrossSectionSearchBarFeatures;
    
    def set_crossSection_barFeatures_isSeaAtWesternSideOfBeach(self, flag):
        self.isSeaAtWesternSideOfBeach = flag;
        
    def get_crossSection_barFeatures_isSeaAtWesternSideOfBeach(self):
        return self.isSeaAtWesternSideOfBeach;
    
    def set_crossSection_barFeatures_relativeOrientation_deg(self, degree):
        self.crossSection_barFeatures_relativeOrientation_deg = degree;
    
    def get_crossSection_barFeatures_relativeOrientation_deg(self):
        return self.crossSection_barFeatures_relativeOrientation_deg;
    
    def set_crossSection_barFeatures_interSpacing_meter(self, meter):
        self.crossSection_barFeatures_interSpacing_meter = meter;
    
    def get_crossSection_barFeatures_interSpacing_meter(self):
        return self.crossSection_barFeatures_interSpacing_meter;
    
    def set_crossSection_barFeatures_min_length_crossSection_meter(self, meter):
        self.crossSection_barFeatures_min_length_crossSection_meter = meter;
    
    def get_crossSection_barFeatures_min_length_crossSection_meter(self):
        return self.crossSection_barFeatures_min_length_crossSection_meter;
    
    def set_crossSection_barFeatures_min_distance_crest_trough_meter(self, meter):
        self.crossSection_barFeatures_min_distance_crest_trough_meter = meter;
    
    def get_crossSection_barFeatures_min_distance_crest_trough_meter(self):
        return self.crossSection_barFeatures_min_distance_crest_trough_meter;
    
    def set_crossSection_barFeatures_min_topo_not_take_into_account_meter(self, meter):
        self.crossSection_barFeatures_min_topo_not_take_into_account_meter = meter;
    
    def get_crossSection_barFeatures_min_topo_not_take_into_account_meter(self):
        return self.crossSection_barFeatures_min_topo_not_take_into_account_meter;
    
    def set_searchIntervalLowReferenceLine(self, tuple_):
        self.searchIntervalLowReferenceLine = tuple_;
    
    def get_searchIntervalLowReferenceLine(self):
        return self.searchIntervalLowReferenceLine;
    
    def set_searchIntervalMeanLowWater(self, tuple_):
        self.searchIntervalMeanLowWater = tuple_;
    
    def get_searchIntervalMeanLowWater(self):
        return self.searchIntervalMeanLowWater;
        
    def set_searchIntervalMeanIntertidal(self, tuple_):
        self.searchIntervalMeanIntertidal = tuple_;
    
    def get_searchIntervalMeanIntertidal(self):
        return self.searchIntervalMeanIntertidal;
    
    def set_searchIntervalMeanHighWater(self, tuple_):
        self.searchIntervalMeanHighWater = tuple_;
    
    def get_searchIntervalMeanHighWater(self):
        return self.searchIntervalMeanHighWater;
    
    def set_searchIntervalHighestAstronomicalTide(self, tuple_):
        self.searchIntervalHighestAstronomicalTide = tuple_;
    
    def get_searchIntervalHighestAstronomicalTide(self):
        return self.searchIntervalHighestAstronomicalTide;
    
    def set_isDoCrossSectionSearchChannelFeatures(self, flag):
        self.isDoCrossSectionSearchChannelFeatures = flag;
    
    def get_isDoCrossSectionSearchChannelFeatures(self):
        return self.isDoCrossSectionSearchChannelFeatures;
    
    def set_crossSection_plotting_frequency(self, freq):
        self.crossSection_plotting_frequency = freq;
    
    def get_crossSection_plotting_frequency(self):
        return self.crossSection_plotting_frequency;
    
    ###########################################################################
    # USER INPUT RELATED TO CHANNEL FEATURES
    ###########################################################################
    def _init_crossSection_channelFeatures(self):
        crossSection_channelFeatures = config.CROSS_SECTION_SEARCH_CHANNEL_FEATURES;
        
        self.set_isDoCrossSectionSearchChannelFeatures(
            crossSection_channelFeatures['DO']
        );
        self.set_crossSection_channelFeatures_relativeOrientation_deg(
            crossSection_channelFeatures['RELATIVE_ORIENTATION_DEGREE']
        );
        self.set_crossSection_channelFeatures_interSpacing_meter(
            crossSection_channelFeatures['INTER-SPACING_METER']
        );
        self.set_crossSection_channelFeatures_min_length_crossSection_meter(
            crossSection_channelFeatures['MINIMUM_LENGTH_CROSSSECTION_METER']
        );
        self.set_crossSection_channelFeatures_minChannelWidth_meter(
            crossSection_channelFeatures['MINIMUM_WIDTH_METER']
        );
        self.set_crossSection_channelFeatures_minChannelDepth_meter(
            crossSection_channelFeatures['MINIMUM_DEPTH_METER']
        );
        self.set_crossSection_channelFeatures_plotting_frequency(
            crossSection_channelFeatures['FREQUENCY_PLOTTING_CROSSSECTIONS']
        );
    
    def set_crossSection_channelFeatures_relativeOrientation_deg(self, degree):
        self.crossSection_channelFeatures_relativeOrientation_deg = degree;
    
    def get_crossSection_channelFeatures_relativeOrientation_deg(self):
        return self.crossSection_channelFeatures_relativeOrientation_deg;
    
    def set_crossSection_channelFeatures_interSpacing_meter(self, meter):
        self.crossSection_channelFeatures_interSpacing_meter = meter;
    
    def get_crossSection_channelFeatures_interSpacing_meter(self):
        return self.crossSection_channelFeatures_interSpacing_meter;
    
    def set_crossSection_channelFeatures_min_length_crossSection_meter(self, meter):
        self.crossSection_channelFeatures_min_length_crossSection_meter = meter;
    
    def get_crossSection_channelFeatures_min_length_crossSection_meter(self):
        return self.crossSection_channelFeatures_min_length_crossSection_meter;
    
    def set_crossSection_channelFeatures_minChannelWidth_meter(self, width):
        self.crossSection_channelFeatures_minChannelWidth = width;
    
    def get_crossSection_channelFeatures_minChannelWidth_meter(self):
        return self.crossSection_channelFeatures_minChannelWidth;
    
    def set_crossSection_channelFeatures_minChannelDepth_meter(self, depth):
        self.crossSection_channelFeatures_minChannelDepth = depth;
    
    def get_crossSection_channelFeatures_minChannelDepth_meter(self):
        return self.crossSection_channelFeatures_minChannelDepth;
    
    def set_crossSection_channelFeatures_plotting_frequency(self, freq):
        self.crossSection_channelFeatures_plotting_frequency = freq;
    
    def get_crossSection_channelFeatures_plotting_frequency(self):
        return self.crossSection_channelFeatures_plotting_frequency;
    
    ###########################################################################
    # USER INPUT RELATED TO POINT FILTERING TO REMOVE ISOLATED POINTS
    ###########################################################################
    def _init_pointFiltering(self):
        point_filtering = config.FILTERING_POINTS;
        self.set_isDoPointFiltering(
            point_filtering['DO']
        );
        self.set_filter_beach_feature_type(
            point_filtering['BEACH_FEATURE_TYPE']
        );
        self.set_filterDimensions_meter(
            point_filtering['FILTER_DIMENSIONS_METER']
        );
        self.set_filterRotation_degree(
            point_filtering['FILTER_ROTATION_DEGREE']
        );
        self.set_filterIterationsOnRaster(
            point_filtering['NUMBER_ITERATIONS_TO_APPLY_FILTER_ON_RASTER']
        );
        self.set_thresholdNeighbourhood_pointCount(
            point_filtering['THRESHOLD_POINT_COUNT_IN_NEIGHBOURHOOD_PCT']
        );
        self.set_thresholdNeighbourhood_topoDiff(
            point_filtering['THRESHOLD_TOPOGRAPHY_DIFFERENCE_IN_NEIGHBOURHOOD']
        );
        self.set_thresholdNeighbourhood_topoDiff_pointCount(
            point_filtering['THRESHOLD_TOPOGRAPHY_DIFFERENCE_IN_NEIGHBOURHOOD_POINT_COUNT_PCT']
        );
    
    def set_isDoPointFiltering(self, flag):
        self.isDoPointFiltering = flag;
    
    def get_isDoPointFiltering(self):
        return self.isDoPointFiltering;
    
    def set_filter_beach_feature_type(self, type_):
        self.filter_beach_feature_type = type_;
    
    def get_filter_beach_feature_type(self):
        return self.filter_beach_feature_type;
    
    def set_filterDimensions_meter(self, dimensions):
        self.filterDimensions_meter = dimensions;
    
    def get_filterDimensions_meter(self):
        return self.filterDimensions_meter;
    
    def set_filterRotation_degree(self, rotation_angle):
        self.filterRotation_degree = rotation_angle;
    
    def get_filterRotation_degree(self):
        return self.filterRotation_degree;
    
    def set_filterIterationsOnRaster(self, number):
        self.filterIterationsOnRaster = number;
    
    def get_filterIterationsOnRaster(self):
        return self.filterIterationsOnRaster;

    def set_thresholdNeighbourhood_pointCount(self, count):
        self.thresholdNeighbourhood_pointCount = count / 100;
    
    def get_thresholdNeighbourhood_pointCount(self):
        return self.thresholdNeighbourhood_pointCount;
    
    def set_thresholdNeighbourhood_topoDiff(self, topography_difference):
        self.thresholdNeighbourhood_topoDiff = topography_difference;
    
    def get_thresholdNeighbourhood_topoDiff(self):
        return self.thresholdNeighbourhood_topoDiff;
    
    def set_thresholdNeighbourhood_topoDiff_pointCount(self, count):
        self.thresholdNeighbourhood_topoDiff_pointCount = count / 100;
    
    def get_thresholdNeighbourhood_topoDiff_pointCount(self):
        return self.thresholdNeighbourhood_topoDiff_pointCount;
    
    ###########################################################################
    # USER INPUT RELATED TO GROUPING POINTS INTO LABELS
    ###########################################################################
    def _init_point_grouping_into_labels(self):
        point_labelling = config.GROUPING_POINTS_IN_LABELS;
        self.set_isDoPointGroupingIntoLabels(
            point_labelling['DO']
        );
        self.set_number_iterations_dilation_m(
            point_labelling['NUMBER_ITERATIONS_DILATION_TO_GROUP_POINTS_METER']
        );
    
    def set_isDoPointGroupingIntoLabels(self, flag):
        self.isDoPointGroupingIntoLabels = flag;
    
    def get_isDoPointGroupingIntoLabels(self):
        return self.isDoPointGroupingIntoLabels;
    
    def set_number_iterations_dilation_m(self, number):
        self.number_iterations_dilation_m = number;
    
    def get_number_iterations_dilation_m(self):
        return self.number_iterations_dilation_m;
    
    ###########################################################################
    # USER INPUT RELATED TO CLUSTERING POINTS
    ###########################################################################
    def _init_point_clustering(self):
        point_clustering = config.GROUPING_POINTS_IN_CLUSTERS;
        self.set_isDoPointClustering(
            point_clustering['DO'])
        self.set_clustering_referenceLine_slope_degree(
            point_clustering['REFERENCE_LINE_SLOPE_DEGREE']
        );
        self.set_clustering_referenceLine_point_coords(
            point_clustering['REFERENCE_LINE_POINT_COORDS']
        );
        self.set_clustering_dbscan_parameters(
            point_clustering['DBSCAN_PARAMETERS']
        );
        self.set_clustering_merge_topographyInterval(
            point_clustering['MERGE_CLUSTERS_TOPOGRAPHY_INTERVAL']
        );
        self.set_clustering_min_spatial_extent_pct(
            point_clustering['MINIMUM_SPATIAL_EXTENT_OF_CLUSTER_PCT']
        );
    
    def set_isDoPointClustering(self, flag):
        self.isDoPointClustering = flag;
    
    def get_isDoPointClustering(self):
        return self.isDoPointClustering;
    
    def set_clustering_referenceLine_slope_degree(self, slope):
        self.clustering_reference_line_slope = slope;
    
    def get_clustering_referenceLine_slope_degree(self):
        return self.clustering_reference_line_slope;
    
    def set_clustering_referenceLine_point_coords(self, point_coords):
        self.clustering_reference_line_point = point_coords;
    
    def get_clustering_referenceLine_point_coords(self):
        return self.clustering_reference_line_point;
    
    def set_clustering_dbscan_parameters(self, params):
        self.clustering_dbscan_eps = params[0];
        self.clustering_dbscan_min_samples = params[1];
    
    def get_clustering_dbscan_eps(self):
        return self.clustering_dbscan_eps;
    
    def get_clustering_dbscan_min_samples(self):
        return self.clustering_dbscan_min_samples;
    
    def set_clustering_merge_topographyInterval(self, interval):
        self.clustering_merge_topography_interval = interval;
    
    def get_clustering_merge_topographyInterval(self):
        return self.clustering_merge_topography_interval;

    def set_clustering_min_spatial_extent_pct(self, percentage):
        self.clustering_min_spatial_extent_pct = percentage;
    
    def get_clustering_min_spatial_extent_pct(self):
        return self.clustering_min_spatial_extent_pct;


class Processing(object):
    def __init__(self, userObject, folderObject, flightObject):
        self.set_userObject(userObject);
        self.set_folderObject(folderObject);
        self.set_flightObject(flightObject);
        
        self.set_rasterObject();
        self.set_targetBounds();
    
    def set_userObject(self, object_):
        self.userObject = object_;
    
    def get_userObject(self):
        return self.userObject;

    def set_flightObject(self, object_):
        self.flightObject = object_;
    
    def get_flightObject(self):
        return self.flightObject;
    
    def set_folderObject(self, object_):
        self.folderObject = object_;
    
    def get_folderObject(self):
        return self.folderObject;

    def set_rasterObject(self, rasterObject=None):
        self.rasterObject = rasterObject;
    
    def get_rasterObject(self):
        return self.rasterObject;
    
    def get_fileExtension(self):
        return '.tif';

    def set_targetBounds(
            self, bounds={'left':-1, 'bottom':-1, 'right':-1, 'top':-1}):
        self.targetBounds = bounds;
    
    def get_targetBounds(self):
        return self.targetBounds;
    
    def set_dtm_res(self, dtm_res):
        self.dtm_res = dtm_res;
    
    def get_dtm_res(self):
        return self.dtm_res;

    def set_topo_mode_res(self, resolution):
        self.topo_mode_res = str(resolution);
    
    def get_topo_mode_res(self):
        return self.topo_mode_res;
    
    def set_flag_smoothing(self, flag):
        self.flag_smoothing = flag;
    
    def get_flag_smoothing(self):
        return self.flag_smoothing;
    
    def set_smooth_res(self, resolution):
        self.smooth_res = str(resolution);
    
    def get_smooth_res(self):
        return self.smooth_res;
        

class Folder(object):     
    def __init__(self):
        self.set_directoryRoot(config.ROOT_DIRECTORY);
    
    def set_flightObject(self, object_):
        self.flightObject = object_;
        
        self.flightName = self.flightObject.get_flightName();
        self.acquisitionPlatform = self.flightObject.get_acquisitionPlatform();
        self.siteName = self.flightObject.get_siteName();
        self.acquisitionDate = self.flightObject.get_acquisitionDate();
        self.flightNumber = self.flightObject.get_flightNumber();
        self.sensorType = self.flightObject.get_sensorType();
    
        self.set_fileName_output_base();
    
        self.set_folderPath_roi();
        self.set_folderPath_inputData();
        self.set_folderPath_processing();
    
        self.set_filePath_shapeFile(
            config.RASTER_PREPROCESSING['CROP_ACCORDING_SHAPE_FILE']
        );
    
    def set_directoryRoot(self, directory):
        self.directoryRoot = directory;
    
    def get_directoryRoot(self):
        return self.directoryRoot;
    
    def set_folderPath_roi(self, directory=None):
        """Defines the directory where to find shapefiles (roi = region of interest)"""
        if directory is None:
            self.folderPath_roi = os.path.join(
                self.directoryRoot,
                'roi',
                self.siteName
            );
        else:
            self.folderPath_roi = directory;
    
    def get_folderPath_roi(self):
        return self.folderPath_roi;
    
    def set_folderPath_inputData(self, directory=None):
        if directory is None:
            self.folderPath_inputData = os.path.join(
                self.directoryRoot,
                self.acquisitionPlatform,
                self.flightName,
                self.sensorType,
                'inputData'
            );
        else:
            self.folderPath_inputData = directory;
    
    def get_folderPath_inputData(self):
        return self.folderPath_inputData;
    
    def get_filePath_inputData(self):
        l_fileNames = [
            x for x in os.listdir(self.folderPath_inputData) if x.endswith('.tif')
        ];
        if len(l_fileNames) != 0:
            filePath = os.path.join(
                self.folderPath_inputData, l_fileNames[0]
            );
            return filePath;
        else:
            logger.warning(
                'Input folder %s does not contain valid geoTIFF file' % self.folderPath_inputData
            );
            raise FileNotFoundError;
    
    def set_filePath_shapeFile(self, file_name=None):
        if file_name is None:
            self.filePath_shapeFile = None;
        else:
            self.filePath_shapeFile = os.path.join(
                self.folderPath_roi, file_name
            );
    
    def get_filePath_shapeFile(self):
        return self.filePath_shapeFile;
    
    def set_fileName_output_base(self, file_name=None):
        if file_name is None:
            self.fileName_output_base = '%s_%s_%s_%s_' % (
                self.acquisitionPlatform,
                self.acquisitionDate,
                self.siteName,
                self.flightNumber
            );
        else:
            self.fileName_output_base = file_name;
    
    def get_fileName_output_base(self):
        return self.fileName_output_base;

    def set_folderPath_processing(self, directory=None):
        if directory is None:
            self.folderPath_processing = os.path.join(
                os.path.split(self.folderPath_inputData)[0],
                'processedData'
            );
        else:
            self.folderPath_processing = directory;
    
    def get_folderPath_processing(self):
        return self.folderPath_processing;

    def get_folderPath_croppedInputData(self, isSmoothed=False):
        if not isSmoothed:
            folderPath = os.path.join(
                self.folderPath_processing, '00_CroppedInputData_NonSmooth'
            );
        else:
            folderPath = os.path.join(
                self.folderPath_processing, '01_CroppedInputData_Smooth'
            ); 
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_crossSections(self):
        folderPath = os.path.join(
            self.folderPath_processing, '10_Cross_Sections'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_topographyMarkers(self):
        folderPath = os.path.join(
            self.folderPath_processing, '11_Topography_Markers'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_relativeExtrema(self):
        folderPath = os.path.join(
            self.folderPath_processing, '12_Crests_and_Troughs'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_relativeExtrema_plot(self):
        folderPath = os.path.join(
            self.get_folderPath_relativeExtrema(),
            'CrossSectionProfiles_Examples'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_barCharacteristics(self):
        folderPath = os.path.join(
            self.folderPath_processing, '13_2D_Bar_Characteristics'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_isolatedRelativeExtremaRemoved(self):
        folderPath = os.path.join(
            self.folderPath_processing, '14_Isolated_Relative_Extrema_Removed'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_fromPointToGroupLabel(self):
        folderPath = os.path.join(
            self.folderPath_processing, '15_From_Point_To_Group_Label'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;
    
    def get_folderPath_fromPointToCluster(self):
        folderPath = os.path.join(
            self.folderPath_processing, '16_From_Point_To_Cluster'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;

    def get_folderPath_timeSeries(self):
        folderPath = os.path.join(
            self.folderPath_processing, '17_Evolution_with_Time'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;

    def get_folderPath_channelLocations(self):
        folderPath = os.path.join(
            self.folderPath_processing, '20_Channel_Locations'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;

    def get_folderPath_channelLocations_plot(self):
        folderPath = os.path.join(
            self.get_folderPath_channelLocations(),
            'CrossSectionProfiles_Examples'
        );
        if not os.path.isdir(folderPath):
            os.makedirs(folderPath);
        return folderPath;

 
class Flight(object):
    def __init__(self):
        self.set_acquisitionPlatform(config.PLATFORM);
        self.set_sensorType(config.SENSOR);
        self.set_siteName(config.SITE_NAME);
    
    def set_flightName(self, name):
        self.flightName = name;
    
        self.set_acquisitionDate();
        self.set_flightNumber();
    
    def get_flightName(self):
        return self.flightName;
    
    def set_acquisitionDate(self, date=None):
        if date is None:
            self.acquisitionDate = self.flightName.split('_')[0];
        else:
            self.acquisitionDate = date; # Format: YYYYMMDD
    
    def get_acquisitionDate(self):
        return self.acquisitionDate;

    def set_flightNumber(self, number=None):
        if number is None:
            self.flightNumber = self.flightName.split('_')[2];
        elif (number >= 0) and (number < 10):
            self.flightNumber = 'F0%s' % number;
        elif (number >=10):
            self.flightNumber = 'F%s' % number;
    
    def get_flightNumber(self):
        return self.flightNumber;

    def set_siteName(self, site):
        self.siteName = site;
    
    def get_siteName(self):
        return self.siteName;

    def set_acquisitionPlatform(self, platform):
        self.acquisitionPlatform = platform;
    
    def get_acquisitionPlatform(self):
        return self.acquisitionPlatform;

    def set_sensorType(self, type_):
        self.sensorType = type_;
    
    def get_sensorType(self):
        return self.sensorType;