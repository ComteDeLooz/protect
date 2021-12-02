# -*- coding: utf-8 -*-

import logging, logging.config
import os
import numpy as np
import pandas as pd

import project
import constants
import raster.platform
import raster.pyraster
import cross_section_analysis.intertidal_bar as cs_bar
import cross_section_analysis.drainage_channel as cs_channel
import cross_section_analysis.water_line as cs_water_line
import cross_section_analysis.utils as cs_utils
import point_analysis.filtering as filtering
import point_analysis.labelling as labelling
import point_analysis.clustering as clustering
import plotting.plot_labelling
import plotting.plot_clustering

# Initialize logging file
logging.config.fileConfig('logging.conf');
logger = logging.getLogger('protect');
logger.info("log protect INFO level");
logger.debug("log protect DEBUG level");


# Generate objects from different container classes in the project.py module 
userObject = project.User();
folderObject = project.Folder();
flightObject = project.Flight();
# Specify the directory where all data related to the mentioned acquisition
#  platform can be found
directoryTarget = os.path.join(
    folderObject.get_directoryRoot(),
    flightObject.get_acquisitionPlatform()
);
logger.info(
    'Generate a list of flights available in the directory under the '\
    + 'acquisition platform folder. Only the flights at the user-specified site '\
    + 'are included'
);
lFlightNames = [
    x for x in os.listdir(directoryTarget)\
    if flightObject.get_siteName() in x
];
for flightName in lFlightNames:
    logger.info(
        '-------------------------------------------------------------------\n'\
        + 'Process data in folder %s' % flightName
    );
    ###########################################################################
    # RASTER PRE-PROCESSING
    ###########################################################################
    flightObject.set_flightName(flightName);
    folderObject.set_flightObject(flightObject);
    if not os.path.isdir(folderObject.get_folderPath_inputData()):
        lFlightNames.remove(flightName);
        logger.warning(
            'Raster file for flight name %s not found based on user-specified information. '\
            + 'Flight name is deleted from list. Processing is continued '\
            + 'with next flight name in the list.' % flightName
        );
        continue;
    # Making object from container class project.Processing 
    processingObject = project.Processing(
        userObject, folderObject, flightObject
    );
    ###########################################################################
    # Read raster data
    logger.info('Read %s data' % flightObject.get_acquisitionPlatform());
    try:
        filePath_inputData = folderObject.get_filePath_inputData();
    except FileNotFoundError:
        logger.warning(
            'Processing is continued with next flight name in the list'
        );
        continue;
    rasterObject = raster.platform.DataAcquiredFromPlatform(
        filePath_inputData,
        folderObject.get_filePath_shapeFile(),
    );
    ###########################################################################
    logger.info('Raster filtering to maintain only the beach area');
    topo_threshold_low = userObject.get_filterBeachAreaOnRaster_topoThresholdLow();
    topo_threshold_high = userObject.get_filterBeachAreaOnRaster_topoThresholdHigh();
    rasterObject.apply_hinterlandFilter(
        threshold_topo=topo_threshold_high
    );
    rasterObject.apply_waterFilter(
        threshold_topo=topo_threshold_low
    );
    rasterObject.apply_beachFilter();    
    ###########################################################################
    if userObject.get_isDoTopoSmoothing():
        logger.info('Smoothing of the topography raster');
        rasterObject.apply_rasterSmooth(
            rasterType='Topography',
            sigma_input=userObject.get_topoSmoothingSigma(),
            sigma_in_pxl=userObject.get_isTopoSmoothingSigmaExpressedInPixels(),
            threshold_topo_low=userObject.get_filterBeachAreaOnRaster_topoThresholdLow(),
            threshold_topo_high=userObject.get_filterBeachAreaOnRaster_topoThresholdHigh()
        );
    ###########################################################################
    # Storing raster object in container object from project.Processing class
    processingObject.set_rasterObject(rasterObject);
    rTransform_window = rasterObject.get_transformWindow(); 
    ###########################################################################
    logger.info('Output the cropped raw (non-smoothed) topography raster');
    # Output the cropped raw (non-smoothed) topography raster
    raster_topo_raw = rasterObject.get_rasterRaw(
        rasterType='Topography', isMasked=False
    );
    filePath_out = os.path.join(
        folderObject.get_folderPath_croppedInputData(isSmoothed=False),
        '%sraw.tif' % folderObject.get_fileName_output_base()
    );
    raster.pyraster.export_geotiff(
        raster_input=raster_topo_raw,
        transform_out=rTransform_window,
        filePath_out=filePath_out,
        dtype_out=np.float32,
        epsg_out=rasterObject.get_epsg(),
        noDataValue=rasterObject.get_noDataValue()
    );
    ###########################################################################
    if userObject.get_isDoTopoSmoothing():
        logger.info('Output the cropped smoothed topography raster');
        raster_topo_smooth = rasterObject.get_rasterSmooth(
            rasterType='Topography', isMasked=False
        );
        filePath_out = os.path.join(
            folderObject.get_folderPath_croppedInputData(isSmoothed=True),
            '%ssmooth.tif' % folderObject.get_fileName_output_base()
        );
        raster.pyraster.export_geotiff(
            raster_input=raster_topo_smooth,
            transform_out=rTransform_window,
            filePath_out=filePath_out,
            dtype_out=np.float32,
            epsg_out=rasterObject.get_epsg(),
            noDataValue=rasterObject.get_noDataValue()
        );
    else:
        logger.warning(
            'User did not activate Gaussian smoothing. '\
            + 'In order to keep workflow operational, the raw raster is copied. '\
            + 'No smoothed raster is exported.'
        );
        raster_topo_smooth = np.copy(raster_topo_raw);
    ###########################################################################
    # FINDING SPECIFIC MARKER POINTS BASED ON TOPOGRAPHY
    # FINDING MAIN INTERTIDAL CRESTS AND TROUGHS
    # CALCULATING 2D INTERTIDAL BAR CHARACTERISTICS
    ###########################################################################    
    absolute_orientation_beach = rasterObject.get_beachOrientation('Degrees');
    logger.info(
        'Absolute orientation of beach is set at %s degree' % str(absolute_orientation_beach)
    );
    if userObject.get_isDoCrossSectionSearchBarFeatures():
        logger.info(
            'Start analysis of cross-sections over the beach area...'
        );
        #######################################################################
        # Get (x,y) coordinates of the start and end points of the
        #   cross-sections with the given orientation
        relative_orientation_crossSection = userObject.get_crossSection_barFeatures_relativeOrientation_deg();
        absolute_orientation_crossSection = np.int(
            absolute_orientation_beach + relative_orientation_crossSection
        );
        # Keep the orientation value in the interval [-90° ; +90°]
        if absolute_orientation_crossSection > 90:
            absolute_orientation_crossSection -= 180;
        elif absolute_orientation_crossSection < -90:
            absolute_orientation_crossSection += 180;
        logger.info('Get the start and end coordinates of the cross-sections');
        crossSections_start_x,\
        crossSections_start_y,\
        crossSections_end_x,\
        crossSections_end_y = cs_utils.get_coordinates_crossSections_refRaster(
            raster_pixel_resolution=rasterObject.get_pixelResolution(),
            raster_shape=rasterObject.get_shape(),
            profile_orientation=absolute_orientation_crossSection,
            profile_spacing_m=userObject.get_crossSection_barFeatures_interSpacing_meter()
        );
        # The function above puts the start of the coordinates on the left/top
        #  of the raster and considers this as the sea part. If this is not the
        #  case (because the sea is at the eastern part of raster), the cross-sections
        #  should be reversed.
        if not userObject.get_crossSection_barFeatures_isSeaAtWesternSideOfBeach():
            logger.info(
                'Reverse cross-section profiles to take into account the sea '\
                'is located at the eastern part of the raster'
            );
            crossSections_start_x_prev = np.copy(crossSections_start_x);
            crossSections_start_y_prev = np.copy(crossSections_start_y);
            crossSections_end_x_prev = np.copy(crossSections_end_x);
            crossSections_end_y_prev = np.copy(crossSections_end_y);
            
            crossSections_start_x = crossSections_end_x_prev;
            crossSections_start_y = crossSections_end_y_prev;
            crossSections_end_x = crossSections_start_x_prev;
            crossSections_end_y = crossSections_start_y_prev;
        #######################################################################
        folderPath_topoMarkers = folderObject.get_folderPath_topographyMarkers();
        filePath_lowReferencePoints = os.path.join(
            folderPath_topoMarkers,
            '%sLowReferencePoints.csv' % folderObject.get_fileName_output_base()
        );
        # If the file path does not exist: do the reference line analysis
        if not os.path.exists(filePath_lowReferencePoints):
            logger.info(
                'Finding water line markers on the raster by analyzing cross-sections'
            );
            cs_waterMarker = cs_water_line.WaterReferenceLinesOnRaster(
                rasterObject.get_pixelResolution(),
                rasterObject.get_coordsRasterOrigin_refEN(),
                raster_topo_raw,
                userObject.get_searchIntervalLowReferenceLine(),
                userObject.get_searchIntervalMeanLowWater(),
                userObject.get_searchIntervalMeanIntertidal(),
                userObject.get_searchIntervalMeanHighWater(),
                userObject.get_searchIntervalHighestAstronomicalTide()
            );
            cs_waterMarker.run_through_crossSections_over_raster(
                crossSections_start_x, crossSections_start_y,
                crossSections_end_x, crossSections_end_y
            );
            pd_lr_raster = cs_waterMarker.get_lowReferencePoints_over_complete_raster();
            if not pd_lr_raster.empty:
                logger.info('Storing points of low reference line');
                pd_lr_raster.to_csv(filePath_lowReferencePoints, sep=';');
            pd_mlw_raster = cs_waterMarker.get_meanLowWaterPoints_over_complete_raster();
            if not pd_mlw_raster.empty:
                logger.info('Storing points of mean low water line');
                pd_mlw_raster.to_csv(os.path.join(
                    folderPath_topoMarkers,
                    '%sMeanLowWaterPoints.csv' % folderObject.get_fileName_output_base()
                ), sep=';');
            pd_mi_raster = cs_waterMarker.get_meanIntertidalPoints_over_complete_raster();
            if not pd_mi_raster.empty:
                logger.info('Storing points of mean intertidal line');
                pd_mi_raster.to_csv(os.path.join(
                    folderPath_topoMarkers,
                    '%sMeanIntertidalPoints.csv' % folderObject.get_fileName_output_base()
                ), sep=';');
            pd_mhw_raster = cs_waterMarker.get_meanHighWaterPoints_over_complete_raster();
            if not pd_mhw_raster.empty:
                logger.info('Storing points of mean high water line');
                pd_mhw_raster.to_csv(os.path.join(
                    folderPath_topoMarkers,
                    '%sMeanHighWaterPoints.csv' % folderObject.get_fileName_output_base()
                ), sep=';');
            pd_hat_raster = cs_waterMarker.get_highestAstronomicalTidePoints_over_complete_raster();
            if not pd_hat_raster.empty:
                logger.info('Storing points of highest astronomical tide line');
                pd_hat_raster.to_csv(os.path.join(
                    folderPath_topoMarkers,
                    '%sHighestAstronomicalTidePoints.csv' % folderObject.get_fileName_output_base()
                ), sep=';');
        else:
            logger.info(
                'File for low reference line exists. '\
                + 'Recalculation of water line markers is skipped.'
            );
        #######################################################################
        # Finding crest-trough markers on the raster
        # Calculating 2D bar characteristics on the cross-section profiles
        folderPath_crestTrough = folderObject.get_folderPath_relativeExtrema();
        filePath_crestPoints = os.path.join(
            folderPath_crestTrough,
            '%sCrestPoints.csv' % folderObject.get_fileName_output_base()
        );
        if not os.path.exists(filePath_crestPoints):
            logger.info(
                'Searching for inflection, crest and trough points on the cross-sections'
            );
            cs_crestTrough = cs_bar.CrestTroughPointsOnRaster(
                raster_pixel_resolution=rasterObject.get_pixelResolution(),
                raster_origin_refEN=rasterObject.get_coordsRasterOrigin_refEN(),
                raster_topo_raw=raster_topo_raw,
                raster_topo_smooth=raster_topo_smooth,
            );
            try:
                pd_lowReferencePoints = pd.read_csv(
                    filePath_lowReferencePoints, sep=';', header=0
                );
            except FileNotFoundError:
                pd_lowReferencePoints = pd.DataFrame();
            # Search for feature points troughout the raster by running
            #  through all cross-section profiles available
            cs_crestTrough.run_through_crossSections_over_raster(
                pd_lowReferencePoints,
                cs_bar.get_minimum_topography_to_take_into_account(userObject),
                userObject.get_crossSection_barFeatures_min_distance_crest_trough_meter(),
                crossSections_start_x, crossSections_start_y,
                crossSections_end_x, crossSections_end_y,
                userObject.get_crossSection_plotting_frequency(),
                folderObject.get_folderPath_relativeExtrema_plot(),
                rasterObject.get_shape(),
                (topo_threshold_low,topo_threshold_high)
            );
            logger.info('Storing all detected inflection points in a CSV-file');
            pd_inflection_raster = cs_crestTrough.get_inflectionPoints_over_complete_raster();
            pd_inflection_raster.to_csv(os.path.join(
                folderPath_crestTrough,
                '%sInflectionPoints.csv' % folderObject.get_fileName_output_base()
            ), sep=';');
            logger.info('Storing all detected crest points in a CSV-file');
            pd_crest_raster = cs_crestTrough.get_crestPoints_over_complete_raster();
            pd_crest_raster.to_csv(filePath_crestPoints, sep=';');
            logger.info('Storing all detected trough points in a CSV-file');
            pd_trough_raster = cs_crestTrough.get_troughPoints_over_complete_raster();
            pd_trough_raster.to_csv(os.path.join(
                folderPath_crestTrough,
                '%sTroughPoints.csv' % folderObject.get_fileName_output_base()
            ), sep=';');
            logger.info('Storing all 2D measurements in a CSV-file');
            pd_bar_characteristics = cs_crestTrough.get_barCharacteristics_over_complete_raster();
            pd_bar_characteristics.to_csv(os.path.join(
                folderObject.get_folderPath_barCharacteristics(),
                '%s2DBarCharacteristics.csv' % folderObject.get_fileName_output_base()
            ), sep=';');
        else:
            logger.info(
                'File for crest points already exists. '\
                + 'Recalculation of inflection, crest and trough points is skipped.')
    
    ###########################################################################
    # SEARCH THE LOCATIONS OF THE DRAINAGE CHANNELS
    ###########################################################################
    if userObject.get_isDoCrossSectionSearchChannelFeatures():
        logger.info('Searching for feature points of drainage channels...');
        relative_orientation_crossSection = userObject.get_crossSection_channelFeatures_relativeOrientation_deg();
        absolute_orientation_crossSection = np.int(
            absolute_orientation_beach + relative_orientation_crossSection
        );
        crossSections_start_x,\
        crossSections_start_y,\
        crossSections_end_x,\
        crossSections_end_y = cs_utils.get_coordinates_crossSections_refRaster(
            raster_pixel_resolution=rasterObject.get_pixelResolution(),
            raster_shape=rasterObject.get_shape(),
            profile_orientation=absolute_orientation_crossSection,
            profile_spacing_m=userObject.get_crossSection_channelFeatures_interSpacing_meter()
        ); 
      
        cs_channelPoints = cs_channel.ChannelPointsOnRaster(
            raster_pixel_resolution=rasterObject.get_pixelResolution(),
            raster_origin_refEN=rasterObject.get_coordsRasterOrigin_refEN(),
            raster_topo_raw=raster_topo_raw,
            raster_topo_smooth=raster_topo_smooth,
            minimum_channel_width_m=userObject.get_crossSection_channelFeatures_minChannelWidth_meter(),
            minimum_channel_depth_m=userObject.get_crossSection_channelFeatures_minChannelDepth_meter()
        );
        cs_channelPoints.run_through_crossSections_over_raster(
            crossSections_start_x, crossSections_start_y,
            crossSections_end_x, crossSections_end_y,
            userObject.get_crossSection_channelFeatures_plotting_frequency(),
            folderObject.get_folderPath_channelLocations_plot(),
            rasterObject.get_shape(),
            (topo_threshold_low,topo_threshold_high)
        );
        _, _,\
        pd_bottom_raster_raw,\
        _, _ = cs_channelPoints.get_channelPoints_over_complete_raster();
        pd_bottom_raster_raw.to_csv(os.path.join(
            folderObject.get_folderPath_channelLocations(),
            '%sChannelBottomPoints_raw.csv' % folderObject.get_fileName_output_base()
        ), sep=';');
    
    ###########################################################################
    # FILTERING OUT ISOLATED RELATIVE EXTREMA POINTS 
    ###########################################################################
    if userObject.get_isDoPointFiltering():
        logger.info('Filtering out isolated relative extrema points');
        # Getting the list of files in the directory containing isolated
        #  relative extrema removed
        check_directory_empty = os.listdir(
            folderObject.get_folderPath_isolatedRelativeExtremaRemoved()
        );
        # Define the type of beach structure in use to filter out points.
        topoFeatureType_toBeFiltered = userObject.get_filter_beach_feature_type();
        # Detecting isolated points is automatically done for the topography feature point
        #  type (e.g. crest or trough) and its strict_condition variant
        pdFilteredPoints = None;
    
        # Walk through all the subdirectories of the selected folder path
        for root, _, files in os.walk(folderObject.get_folderPath_relativeExtrema()):
            # If the directory referring to the isolated relative extrema is not empty,
            #  no recalculation of the isolated relative points will be done 
            if len(check_directory_empty) != 0:
                logger.warning(
                    'Directory storing the removed isolated relative extrema '\
                    + 'is not empty, no recalculation of isolated points will be done')
                break;
            # Only consider files in the main directory of relative extrema results
            if root != folderObject.get_folderPath_relativeExtrema():
                continue;
            # Step 1: filter out isolated points of the selected topography feature type 
            logger.info(
                'Working on points of the type %s' % topoFeatureType_toBeFiltered
            );
            for file in files:
                # Get beach structure type out of the name of the feature point CSV
                # file name convention: e.g. LiDAR_date_location_Fxx_CrestPoints.csv
                topoFeatureType = file.split('_')[4][:-4];
                if topoFeatureType != topoFeatureType_toBeFiltered:
                    continue;
                filePath_pointCloud = os.path.join(
                    folderObject.get_folderPath_relativeExtrema(), file
                );
                pc_nf = filtering.NeighbourhoodFilter(
                    filePath_pointCloud,
                    rasterObject.get_pixelResolution(),
                    rasterObject.get_shape()
                );
                filter_rotationAngle_deg = userObject.get_filterRotation_degree();
                # If no rotation angle of the filter window is specified,
                #  select the main orientation of the beach area
                if filter_rotationAngle_deg is None:
                    filter_rotationAngle_deg = absolute_orientation_beach;
                pc_nf.set_filter_characteristics(
                    neighbourhoodDimensions_m=userObject.get_filterDimensions_meter(),
                    rotationAngle_deg=absolute_orientation_beach,
                    number_iterations=userObject.get_filterIterationsOnRaster(),
                    threshold_point_count=userObject.get_thresholdNeighbourhood_pointCount(),
                    threshold_topo_difference=userObject.get_thresholdNeighbourhood_topoDiff(),
                    threshold_topo_difference_point_count=userObject.get_thresholdNeighbourhood_topoDiff_pointCount()
                );
                pc_nf.apply_rectangular_filter();
                
                # Get and store the raster with points that were not
                #  filtered out
                pdFilteredPoints_base = pc_nf.get_filteredPointsInTableForm();
                pdFilteredPoints_base.to_csv(os.path.join(
                    folderObject.get_folderPath_isolatedRelativeExtremaRemoved(),
                    '%s_filtered.csv' % file[:-4]
                ), sep=';');
                # Get and store the information about the applied filtering
                pdFilteredPointsInfo = pc_nf.get_filterInfoTable();
                pdFilteredPointsInfo.to_csv(os.path.join(
                    folderObject.get_folderPath_isolatedRelativeExtremaRemoved(),
                    '%s_filteredInfo.csv' % file[:-4]
                ), sep=';');

            # Step 2: filter out point ID's for the other structure types
            #  (e.g. inflection & trough) based on what is filtered out on
            #  the target beach structure type (e.g. crest) 
            for file in files:
                # Get beach structure type out of the name of the feature point CSV
                # file name convention: e.g. LiDAR_date_location_Fxx_CrestPoints.csv
                beachFeatureType = file.split('_')[4][:-4];
                if beachFeatureType == topoFeatureType_toBeFiltered:
                    continue;
                logger.info(
                    'Working on points of the type %s' % beachFeatureType
                );
                # Get a list of all point ID's which should be maintained (not filtered out)
                npFilteredPoints_pointID = pdFilteredPoints_base.loc[
                    :,'Point_ID'
                ].values;
                # Load the file based on the data that should be filtered.
                filePath_pointCloud = os.path.join(
                    folderObject.get_folderPath_relativeExtrema(), file
                );
                pdPointCloud = pd.read_csv(
                    filePath_pointCloud, sep=';', index_col=0
                );
                pdFilteredPoints = pdPointCloud.loc[
                    pdPointCloud['Point_ID'].isin(npFilteredPoints_pointID), :
                ];
                pdFilteredPoints.to_csv(os.path.join(
                    folderObject.get_folderPath_isolatedRelativeExtremaRemoved(),
                    '%s_filtered.csv' % file[:-4]
                ), sep=';');
    
    ###########################################################################
    # GROUP FEATURE POINTS INTO COMMON LABELS
    ###########################################################################
    if userObject.get_isDoPointGroupingIntoLabels():
        logger.info('Start grouping of intertidal bar feature points into individual bars...');
        folderPath_input = folderObject.get_folderPath_isolatedRelativeExtremaRemoved();
        fileName_extension = '_filtered.csv';
        if len(os.listdir(folderPath_input)) == 0:
            os.rmdir(folderPath_input);
            folderPath_input = folderObject.get_folderPath_relativeExtrema();
            fileName_extension = '.csv'
        logger.info('Use input files in directory %s' % folderPath_input);
        folderPath_output = folderObject.get_folderPath_fromPointToGroupLabel();
        npFeaturePointNames = np.array([
            'CrestPoints', 'TroughPoints', 'InflectionPoints',
        ]);
        pixel_resolution = rasterObject.get_pixelResolution();
        pixel_resolution_median = np.median(
            np.absolute(pixel_resolution)
        );
        l_binaryDilation_numberIterations = np.asarray(
            userObject.get_number_iterations_dilation_m() / pixel_resolution_median,
            dtype=np.int
        );
        for npIx, feature_point_name in np.ndenumerate(npFeaturePointNames):
            filePath_pointCloud_csv = os.path.join(folderPath_input, '%s%s%s' % (
                folderObject.get_fileName_output_base(),
                feature_point_name,
                fileName_extension
            ));
            pdPointCloud = pd.read_csv(
                filePath_pointCloud_csv, sep=';', index_col=0
            );
            # If the raster contains less than x feature points, skip the analysis
            if pdPointCloud.shape[0] < constants.MINIMUM_NUMBER_POINTS_ON_RASTER_TO_DO_ANALYSIS:
                continue;
            cptib = labelling.FromPointToGroupLabel(
                filePath_pointCloud_csv,
                rasterObject.get_shape(),
                pixel_resolution_median
            );
            cptib.group_points_in_labels(
                number_iterations=l_binaryDilation_numberIterations[npIx]
            );
            cptib.calc_intertidalBarCharacteristics();
            pdPointCloud = cptib.get_appendedInformation();
            fileName = '%s%s_label.csv' % (
                folderObject.get_fileName_output_base(),
                feature_point_name
            );
            pdPointCloud.to_csv(
                os.path.join(folderPath_output, fileName), sep=';'
            );  
            plotting.plot_labelling.plot_fromPointToGroupLabel(
                raster_data=rasterObject.get_rasterSmooth('Topography', True),
                raster_shape=rasterObject.get_shape(),
                raster_pixel_resolution=rasterObject.get_pixelResolution(),
                raster_origin_refEN=rasterObject.get_coordsRasterOrigin_refEN(),
                filePath_pdPointCloud=os.path.join(folderPath_output, fileName),
                topo_threshold_low=userObject.get_filterBeachAreaOnRaster_topoThresholdLow(),
                topo_threshold_high=userObject.get_filterBeachAreaOnRaster_topoThresholdHigh()
            );
    
    ###########################################################################
    # CLUSTERING OF POINT CLOUD
    ###########################################################################
    if userObject.get_isDoPointClustering():
        logger.info('Start clustering of intertidal bar feature points...');
        folderPath_input = folderObject.get_folderPath_fromPointToGroupLabel();
        fileName_extension = '_label.csv';
        if len(os.listdir(folderPath_input)) == 0:
            os.rmdir(folderPath_input);
            folderPath_input = folderObject.get_folderPath_isolatedRelativeExtremaRemoved();
            fileName_extension = '_filtered.csv';
        if len(os.listdir(folderPath_input)) == 0:
            os.rmdir(folderPath_input);
            folderPath_input = folderObject.get_folderPath_relativeExtrema();
            fileName_extension = '.csv'; 
        logger.info('Use input files in directory %s' % folderPath_input);
        folderPath_output = folderObject.get_folderPath_fromPointToCluster();
        npFeaturePointNames = np.array([
            'CrestPoints', 'TroughPoints', 'InflectionPoints',
        ]);
        for npIx, feature_point_name in np.ndenumerate(npFeaturePointNames):
            logger.info('Clustering %s' % feature_point_name);
            filePath_pointCloud_csv = os.path.join(folderPath_input, '%s%s%s' % (
                folderObject.get_fileName_output_base(),
                feature_point_name,
                fileName_extension
            ));
            try:
                pdPointCloud = pd.read_csv(
                    filePath_pointCloud_csv, sep=';', index_col=0
                );
            except FileNotFoundError:
                logger.warning(
                    'File with %s could not be found '\
                    '- continue with next intertidal bar feature points' % feature_point_name
                );
                continue;
            # If the raster contains less than x feature points, skip the analysis
            if pdPointCloud.shape[0] < constants.MINIMUM_NUMBER_POINTS_ON_RASTER_TO_DO_ANALYSIS:
                logger.warning(
                    'The number of intertidal bar feature points on the raster '\
                    + 'is insufficient to do analysis'
                );
                continue;
            # In PROTECT script, the orientation is defined (+) in clockwise direction
            reference_line_slope_degree = userObject.get_clustering_referenceLine_slope_degree();
            if reference_line_slope_degree is None:
                reference_line_slope_degree = -absolute_orientation_beach;
            clusterObject = clustering.FromPointToCluster(
                pdPointCloud,
                reference_line_slope=np.tan(np.radians(reference_line_slope_degree)),
                reference_line_point=userObject.get_clustering_referenceLine_point_coords(),
                merge_clusters_interval_width_m=userObject.get_clustering_merge_topographyInterval()
            );
            clusterObject.apply_clustering_dbscan(
                userObject.get_clustering_dbscan_eps(),
                userObject.get_clustering_dbscan_min_samples()
            );
            logger.info('Merge clusters with similar topography values')
            clusterObject.merge_clusters_with_comparable_topography(
                userObject.get_clustering_merge_topographyInterval()
            );
            logger.info(
                'Reorganizing clusters found to make sure only relevant clusters '\
                + 'are maintained and having an ordering number similar to other acquisition dates')
            length_beach_area_pxl = rasterObject.get_beachEstimatedLength();
            spacing_between_profiles_meter = userObject.get_crossSection_barFeatures_interSpacing_meter();
            pixel_resolution_x, pixel_resolution_y = rasterObject.get_pixelResolution();
            pixel_resolution_mean = np.mean(np.absolute([
                pixel_resolution_x, pixel_resolution_y
            ]));
            spacing_between_profiles_pxl = spacing_between_profiles_meter\
                                           / pixel_resolution_mean; 
            clusterObject.reorganize_clusters_found(
                length_beach_area_pxl,
                spacing_between_profiles_pxl,
                min_spatial_extent_pct=userObject.get_clustering_min_spatial_extent_pct()
            );
            pdPointCloud_update = clusterObject.get_point_cloud_with_appended_information();
            logger.info(
                'Plotting clustering results of %s' % feature_point_name
            );
            fileName = '%s%s_cluster.csv' % (
                folderObject.get_fileName_output_base(),
                feature_point_name
            );
            pdPointCloud_update.to_csv(
                os.path.join(folderPath_output, fileName), sep=';'
            );  
            plotting.plot_clustering.plot_fromPointToCluster(
                raster_data=rasterObject.get_rasterSmooth('Topography', True),
                raster_shape=rasterObject.get_shape(),
                raster_pixel_resolution=rasterObject.get_pixelResolution(),
                raster_origin_refEN=rasterObject.get_coordsRasterOrigin_refEN(),
                filePath_pdPointCloud=os.path.join(folderPath_output, fileName),
                topo_threshold_low=userObject.get_filterBeachAreaOnRaster_topoThresholdLow(),
                topo_threshold_high=userObject.get_filterBeachAreaOnRaster_topoThresholdHigh(),
                reference_line_slope=np.tan(np.radians(reference_line_slope_degree)),
                reference_line_point=userObject.get_clustering_referenceLine_point_coords()
            );
            filePath_barChars = os.path.join(
                folderObject.get_folderPath_barCharacteristics(),
                '%s2DBarCharacteristics.csv' % folderObject.get_fileName_output_base()
            );
            try:
                pdBarChars = pd.read_csv(filePath_barChars, sep=';');
            except FileNotFoundError:
                logger.warning(
                    'File with intertidal bar measurements could not be found. '\
                    + 'No merge possible with clustering file to do time series analysis'
                );
                continue;
            pdPointCloud_merge = pdPointCloud_update.merge(
                pdBarChars, how='outer',
                left_on='Point_ID', right_on='Crest_PointID'
            );
            fileName = '%s%s_merge.csv' % (
                folderObject.get_fileName_output_base(),
                feature_point_name
            );
            pdPointCloud_merge.to_csv(os.path.join(
                folderObject.get_folderPath_timeSeries(),
                fileName
            ), sep=';');
