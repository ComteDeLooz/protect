# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy import signal

import config, constants
import cross_section_analysis.utils as cs_utils
import plotting.plot_intertidal_bar as plt_ib


def get_minimum_topography_to_take_into_account(userObject):
    threshold_min_topo_valid = userObject.get_crossSection_barFeatures_min_topo_not_take_into_account_meter();
    if threshold_min_topo_valid is None:
        threshold_min_topo_valid = userObject.get_filterBeachAreaOnRaster_topoThresholdLow();
    if threshold_min_topo_valid is None:
        threshold_min_topo_valid = userObject.get_searchIntervalLowReferenceLine();
    if threshold_min_topo_valid is None:
        threshold_min_topo_valid = -5; # Default value
    return threshold_min_topo_valid;


class CrestTroughPointsOnRaster():
    """
    Processing of cross-section profiles to find maximum of bars (crest)
     and minimum of troughs on the input raster by analyzing cross-sections.
    Measurements of the intertidal bars along the cross-section profiles (2D)
     is also performed.
    """
    def __init__(
            self,
            raster_pixel_resolution,
            raster_origin_refEN,
            raster_topo_raw,
            raster_topo_smooth):
        """
        Parameters
        ----------
        raster_pixel_resolution (tuple of float)
            (pixel_resolution_x, pixel_resolution_y)
        raster_origin_refEN (tuple of float)
            Coordinates (E,N) of the origin of the image (top left)
        raster_topo_raw (numpy.ndarray of float)
            2D raster data of raw (non-smoothed) topography
        raster_topo_smooth (numpy.ndarray of float)
            2D raster data of smoothed topography
        """
        self.rPixelResolution_x,\
        self.rPixelResolution_y = raster_pixel_resolution;
        self.rPixelResolution_mean = np.mean(
            np.absolute([self.rPixelResolution_x, self.rPixelResolution_y])
        );
        self.rEasting_origin, self.rNorthing_origin = raster_origin_refEN;
        self.raster_topo_raw = raster_topo_raw;
        self.raster_topo_smooth = raster_topo_smooth;
        
        self.init_featurePoints_on_raster();

    def _isCrossSectionLengthEnough(self, crossSectionProfile_topo):
        """
        If a cross-section is too short, no relevant data are obtained.
            Threshold: 2 x maximum of the minimum distance between subsequent
            extrema as specified by the user in the config.py file.
        """
        self.isCrossSectionLengthEnough = True;
        crossSectionLength_threshold_pxl = config.CROSS_SECTION_SEARCH_BAR_FEATURES['MINIMUM_LENGTH_CROSSSECTION_METER']\
                                           / self.rPixelResolution_mean;
        if crossSectionProfile_topo[crossSectionProfile_topo > 0].size \
                <= crossSectionLength_threshold_pxl:
            self.isCrossSectionLengthEnough = False;  
  
    def _init_featurePoints_on_crossSectionProfile(self):
        """
        Initialisation of new arrays that contain features found over the selected cross-section.
         A feature point can be an inflection, crest or trough point.
        Re-initializes for each new cross-section
        """
        self.np_inflection_ixCrossSection = np.array([], dtype=np.int);
        self.np_inflection_pointID = np.array([], dtype=np.int);
        
        self.np_crest_ixCrossSection = np.array([], dtype=np.int);
        self.np_crest_topo_raw = np.array([], dtype=np.float);
        
        self.np_trough_ixCrossSection = np.array([], dtype=np.int);
        self.np_trough_topo_raw = np.array([], dtype=np.float);
        
        self.np_crestTrough_pointID = np.array([], dtype=np.int);
        
    def init_featurePoints_on_raster(self):
        """
        Initialisation of all relevant numpy arrays as empty arrays.
         These arrays store all information regarding position and topography
         of the feature points (inflection, crest or trough points)
        """
        #######################################################################
        # Point information for the different feature points found
        # Unique identification number over the complete raster
        # Crest & trough points linked to their inflection point all have the
        #  same ID number
        # It is possible inflection points exist where no crest / trough point
        #  is allocated to (see further in class)
        self.featurePoint_id = 100000;
        # Storing the ID numbers for each inflection point
        self.np_raster_inflection_ID = np.array([], dtype=np.int);
        # Storing the ID numbers for each crest / trough point
        self.np_raster_crestTrough_ID = np.array([], dtype=np.int);
        
        # x-coordinate of crest point positions with reference to the raster (x = columns)
        self.np_raster_crest_x = np.array([], dtype=np.int);
        # y-coordinate of crest point positions with reference to the raster (y = rows)
        self.np_raster_crest_y = np.array([], dtype=np.int);
        # x-coordinate of trough point positions with reference to the raster (x = columns)
        self.np_raster_trough_x = np.array([], dtype=np.int);
        # y-coordinate of trough point positions with reference to the raster (y = rows)
        self.np_raster_trough_y = np.array([], dtype=np.int);
        # x-coordinate of inflection point positions with reference to the raster (x = columns)
        self.np_raster_inflection_x = np.array([], dtype=np.int);
        # y-coordinate of inflection point positions with reference to the raster (y = rows)
        self.np_raster_inflection_y = np.array([], dtype=np.int);
        
        # Coordinates of point positions with reference to the cross-section profile
        self.np_raster_crest_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_trough_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_inflection_ixCrossSection = np.array([], dtype=np.int);
        
        # Topography value of the non-smoothed (raw) raster on the selected point positions
        self.np_raster_crest_topo = np.array([], dtype=np.float32);
        self.np_raster_trough_topo = np.array([], dtype=np.float32);
        self.np_raster_inflection_topo = np.array([], dtype=np.float32);
        
        #######################################################################
        # Cross-section information per inflection point found
        # Unique identification number of the cross-section
        self.np_raster_inflection_crossSection_ID = np.array([], dtype=np.int);
        # Starting coordinate (x-axis) of the cross-section
        self.np_raster_inflection_crossSection_start_x = np.array([], dtype=np.int);
        # Starting coordinate (y-axis) of the cross-section
        self.np_raster_inflection_crossSection_start_y = np.array([], dtype=np.int);
        # Ending coordinate (x-axis) of the cross-section
        self.np_raster_inflection_crossSection_end_x = np.array([], dtype=np.int);
        # Ending coordinate (y-axis) of the cross-section
        self.np_raster_inflection_crossSection_end_y = np.array([], dtype=np.int);
        
        #######################################################################
        # Cross-section information per crest/trough point found
        self.np_raster_crestTrough_crossSection_ID = np.array([], dtype=np.int);
        
        #######################################################################
        # 2D array to store all characteristics of the intertidal bars
        #  One row in the array corresponds to one intertidal bar linked to
        #  one crest point
        self.np_raster_barCharacteristics = np.array([[]]);

    def _find_inflectionPoints_on_crossSectionProfile(self):
        """
        Searching for inflection points on the smoothed topography cross-section profile.
            * Second derivative of the smoothed topography cross-section
            * Median filtering
            * Look for the location in the cross-section where the sign of the
                second derivative changes sign
            * Maintain only the locations where the sign changes from - to +
            * Sort the numpy array from low to high (lower to higher topography)
        """
        # 1. Calculate the second derivative of the smoothed cross-section profile
        profile_topo_firstDeriv = np.gradient(
            self.crossSectionProfile_topo_smooth
        );
        profile_topo_secondDeriv = np.gradient(
            profile_topo_firstDeriv
        );
        
        # 2. Median filtering the second derivative of the topography profile
        medianFilter_kernelSize_pxl = np.int(np.around(
            constants.BAR['INFLECTION_POINTS_MEDIAN_FILTER_KERNEL_METER']\
            / self.rPixelResolution_mean
        ));
        # Make sure the filter kernel is an uneven number
        if np.mod(medianFilter_kernelSize_pxl, 2) == 0:
            medianFilter_kernelSize_pxl += 1;
        profile_topo_secondDeriv_smooth = signal.medfilt(
            profile_topo_secondDeriv, kernel_size=medianFilter_kernelSize_pxl
        );
        
        # 3. Search for positions on the smoothed second derivative of the topography profile
        #   where the curve goes from negative to positive
        # The idea is to only find inflection points on the landward side of the
        #   intertidal bar  
        np_inflection_ixCrossSection = np.where(
            np.diff(np.sign(profile_topo_secondDeriv_smooth)) > 0
        )[0];
        # Sort the locations (indices) of the points found on the cross-section profile
        self.np_inflection_ixCrossSection = np.sort(
            np_inflection_ixCrossSection
        );

    def _append_inflection_from_crossSection_to_raster(self):
        """
        Store the inflection points found on the cross-section into the large
         array for the complete raster
        """
        np_inflection_x = self.cc[
            self.np_inflection_ixCrossSection
        ];
        np_inflection_y = self.rr[
            self.np_inflection_ixCrossSection
        ];
        np_inflection_topo = self.raster_topo_raw[
            np_inflection_y, np_inflection_x
        ];
        
        self.np_raster_inflection_ID = np.append(
            self.np_raster_inflection_ID, self.np_inflection_pointID
        );
        self.np_raster_inflection_x = np.append(
            self.np_raster_inflection_x, np_inflection_x
        );
        self.np_raster_inflection_y = np.append(
            self.np_raster_inflection_y, np_inflection_y
        );
        self.np_raster_inflection_topo = np.append(
            self.np_raster_inflection_topo, np_inflection_topo
        );
        self.np_raster_inflection_ixCrossSection = np.append(
            self.np_raster_inflection_ixCrossSection,
            self.np_inflection_ixCrossSection
        );
            
        numberOfInflectionPoints = self.np_inflection_ixCrossSection.size;
        self.np_raster_inflection_crossSection_ID = np.append(
            self.np_raster_inflection_crossSection_ID,
            np.full((numberOfInflectionPoints,), self.raster_crossSection_count, dtype=np.int)
        );


    def _find_crestOrTroughPoints_on_crossSectionProfile(
            self,
            index_on_crossSectionProfile_thresholdLower,
            threshold_min_topo_valid,
            min_distance_crest_trough_m):
        """Searches the positions of crest and trough points on the cross-section profile"""
        # Loop over the positions of the inflection points found for the specified
        #   cross-section
        for (inflection_pointCount,), inflection_ixCrossSection in np.ndenumerate(
                self.np_inflection_ixCrossSection):
            # Add for the given inflection point the unique feature point ID
            self.np_inflection_pointID = np.append(
                self.np_inflection_pointID, self.featurePoint_id
            );
            self.featurePoint_id += 1;
            
            # If the position of the inflection point is closer to the sea
            #   than the low reference point, no according crest and trough points
            #   are searched for
            if inflection_ixCrossSection <= index_on_crossSectionProfile_thresholdLower:
                continue;
            # Determine the end of the interval where useful data are available
            #   to find crest and trough points.
            # Option 1: the end is set where the smoothed topography from the
            #   inflection point onwards is highest
            end_ixCrossSection = inflection_ixCrossSection + np.argmax(
                self.crossSectionProfile_topo_smooth[inflection_ixCrossSection:]
            );
            # Option 2: the end is set where the smoothed topography from the
            #   inflection point onwards is the first time negative  
            npNegativeValue_ixCrossSection = np.argwhere(
                self.crossSectionProfile_topo_smooth[inflection_ixCrossSection:]\
                < threshold_min_topo_valid
            );
            if npNegativeValue_ixCrossSection.size != 0:
                negativeValue_ixCrossSection = inflection_ixCrossSection + np.amin(
                    npNegativeValue_ixCrossSection
                );
            else:
                negativeValue_ixCrossSection = end_ixCrossSection + 1;
            # If the index of the inflection point is equal or larger than
            #   the index of the highest point in the topography in the cross section
            #   or than the index of the first no data value, do not consider this in the analysis.
            # The inflection point is considered outside the target area on the beach
            if inflection_ixCrossSection >= end_ixCrossSection\
                    or inflection_ixCrossSection >= negativeValue_ixCrossSection:
                continue;
            ###################################################################
            # Search for crest points
            crest_searchInterval = self.crossSectionProfile_topo_raw[
                index_on_crossSectionProfile_thresholdLower\
                :inflection_ixCrossSection
            ];
            try:
                crest_relativePosition = np.argmax(crest_searchInterval);
            except ValueError:
                continue;
            crest_ixCrossSection = index_on_crossSectionProfile_thresholdLower\
                                   + crest_relativePosition;
            # Topography of the crest point is based on raw (non-smoothed) topography
            crest_topo_raw = self.crossSectionProfile_topo_raw[
                crest_ixCrossSection
            ];
            ###################################################################
            # Search for trough points
            endTroughSearch_ixCrossSection = np.amin([
                end_ixCrossSection, negativeValue_ixCrossSection
            ]);
            trough_searchInterval = self.crossSectionProfile_topo_raw[
                inflection_ixCrossSection:endTroughSearch_ixCrossSection
            ];
            try:
                trough_relativePosition = np.argmin(trough_searchInterval);
            except ValueError:
                continue;
            trough_ixCrossSection = inflection_ixCrossSection\
                                    + trough_relativePosition;
            # Topography of the trough point is based on raw (non-smoothed) topography
            trough_topo_raw = self.crossSectionProfile_topo_raw[
                trough_ixCrossSection
            ];
            ###################################################################
            # Checking whether found crest and trough points fulfill some basic conditions
            
            # For the same inflection point, the crest point should be higher
            #   in topography compared to the trough point
            if crest_topo_raw <= trough_topo_raw:
                continue;
            # The positions of the crest and trough point belonging to the same
            #   inflection point should be at least 3 meters away from each other
            #   This to avoid that crest and trough points fall together in the
            #   same point
            dist_crest_trough = trough_ixCrossSection - crest_ixCrossSection;
            dist_threshold = min_distance_crest_trough_m / self.rPixelResolution_mean;
            if dist_crest_trough <= dist_threshold:
                continue;
            
            crest_ixCrossSection = index_on_crossSectionProfile_thresholdLower\
                                   + crest_relativePosition;
            self.np_crest_ixCrossSection = np.append(
                self.np_crest_ixCrossSection, crest_ixCrossSection
            );
            self.np_crest_topo_raw = np.append(
                self.np_crest_topo_raw, crest_topo_raw
            );
            self.np_trough_ixCrossSection = np.append(
                self.np_trough_ixCrossSection, trough_ixCrossSection
            );
            self.np_trough_topo_raw = np.append(
                self.np_trough_topo_raw, trough_topo_raw
            );
            # The crest and trough points together with the corresponding inflection
            #   points all have the same feature point ID.
            self.np_crestTrough_pointID = np.append(
                self.np_crestTrough_pointID, self.featurePoint_id-1
            );

    def _remove_duplicates_localExtremaPoints(self, isDuplicatesCrestPoints=True):
        """
        Decimates the crest AND trough point arrays based on duplicates found for
            OR crest points OR trough points
        
        Parameters
        ----------
        isDuplicatesCrestPoints (bool)
            Selects the type of feature point to use to check for duplicates
                If True: use crest points (default)
                If False: use trough points
        """
        if isDuplicatesCrestPoints:
            np_duplicateArray_ixCrossSection = self.np_crest_ixCrossSection;
        else:
            np_duplicateArray_ixCrossSection = self.np_trough_ixCrossSection;
        
        # new array with unique elements,\
        # Indices in the old array of the (unique) elements of the new array,\
        # Indices in the new array of all elements of the old array,\
        # Number of times the (unique) elements in the new array are counted in the old array
        np_extrema_ixCrossSection_unique,\
        np_extrema_indices_unique,\
        np_extrema_inverse_unique,\
        np_extrema_counts_unique = np.unique(
            np_duplicateArray_ixCrossSection,
            return_index=True, return_inverse=True, return_counts=True);
        
        np_uniqueArray_ixCrossSection_ixOriginalArray = np.array(
            [], dtype=np.int
        );
        for ix_uniqueArray, counts in np.ndenumerate(np_extrema_counts_unique):            
            # Find all the positions (indices) in the original array that correspond
            #   to the selected element in the unique array
            extrema_ixCrossSection_ixOriginalArray = np.argwhere(
                np_extrema_inverse_unique == ix_uniqueArray
            )[:,0];
            
            # In case there was only one crest/trough point found at the position
            if extrema_ixCrossSection_ixOriginalArray.size == 1:
                np_uniqueArray_ixCrossSection_ixOriginalArray = np.append(
                    np_uniqueArray_ixCrossSection_ixOriginalArray,
                    extrema_ixCrossSection_ixOriginalArray
                );
                continue;
            # Find duplicates based on the crest points
            if isDuplicatesCrestPoints:
                ixDecisionPoint = np.argmin(
                    self.np_trough_topo_raw[extrema_ixCrossSection_ixOriginalArray]
                );
            # Find duplicates based on the trough points
            else:
                ixDecisionPoint = np.argmax(
                    self.np_crest_topo_raw[extrema_ixCrossSection_ixOriginalArray]
                );
            np_uniqueArray_ixCrossSection_ixOriginalArray = np.append(
                np_uniqueArray_ixCrossSection_ixOriginalArray,
                extrema_ixCrossSection_ixOriginalArray[ixDecisionPoint]
            );
        
        # Decimate the arrays by removing duplicates
        self.np_crest_ixCrossSection = self.np_crest_ixCrossSection[
            np_uniqueArray_ixCrossSection_ixOriginalArray
        ];
        self.np_crest_topo_raw = self.np_crest_topo_raw[
            np_uniqueArray_ixCrossSection_ixOriginalArray
        ];
        self.np_trough_ixCrossSection = self.np_trough_ixCrossSection[
            np_uniqueArray_ixCrossSection_ixOriginalArray
        ];
        self.np_trough_topo_raw = self.np_trough_topo_raw[
            np_uniqueArray_ixCrossSection_ixOriginalArray
        ];
        self.np_crestTrough_pointID = self.np_crestTrough_pointID[
            np_uniqueArray_ixCrossSection_ixOriginalArray
        ];     

    def _append_crestTrough_from_crossSection_to_raster(self):       
        """
        Store the crest/trough points found on the cross-section into the large
         array for the complete raster
        """
        self.np_raster_crestTrough_ID = np.append(
            self.np_raster_crestTrough_ID, self.np_crestTrough_pointID
        );

        self.np_crest_x = self.cc[
            self.np_crest_ixCrossSection
        ];
        self.np_crest_y = self.rr[
            self.np_crest_ixCrossSection
        ];
        np_crest_topo = self.raster_topo_raw[
            self.np_crest_y, self.np_crest_x
        ];         
        self.np_raster_crest_x = np.append(
            self.np_raster_crest_x, self.np_crest_x
        );
        self.np_raster_crest_y = np.append(
            self.np_raster_crest_y, self.np_crest_y
        );
        self.np_raster_crest_topo = np.append(
            self.np_raster_crest_topo, np_crest_topo
        );
        self.np_raster_crest_ixCrossSection = np.append(
            self.np_raster_crest_ixCrossSection, self.np_crest_ixCrossSection
        );
        
        self.np_trough_x = self.cc[
            self.np_trough_ixCrossSection
        ];
        self.np_trough_y = self.rr[
            self.np_trough_ixCrossSection
        ];
        np_trough_topo = self.raster_topo_raw[
            self.np_trough_y, self.np_trough_x
        ];       
        self.np_raster_trough_x = np.append(
            self.np_raster_trough_x, self.np_trough_x
        );
        self.np_raster_trough_y = np.append(
            self.np_raster_trough_y, self.np_trough_y
        );
        self.np_raster_trough_topo = np.append(
            self.np_raster_trough_topo, np_trough_topo
        );
        self.np_raster_trough_ixCrossSection = np.append(
            self.np_raster_trough_ixCrossSection, self.np_trough_ixCrossSection);
        
        numberOfCrestTroughPoints = self.np_crest_ixCrossSection.size;
        self.np_raster_crestTrough_crossSection_ID = np.append(
            self.np_raster_crestTrough_crossSection_ID,
            np.full((numberOfCrestTroughPoints,), self.raster_crossSection_count, dtype=np.int)
        );

    def _barCharacteristics_on_crossSectionProfile(self):
        """
        Based on the found crest and trough points on the considered cross-section,
         a bar characteristics analysis is performed.
        Measurements take into account ALL points found on the cross-section.
        """
        # Loop through all crest point positions found for the considered cross-section
        #  (sorted from lowest to highest topography)
        for (pointCount,), crest_ixCrossSection in np.ndenumerate(
                self.np_crest_ixCrossSection):
            crest_ID = self.np_crestTrough_pointID[pointCount];
            crest_topo = self.crossSectionProfile_topo_raw[
                crest_ixCrossSection
            ];
            # Initialize all parameters to be measured
            trough_ID_seaward = -1;
            barWidth_total = -1;
            barSymmetry = -1;
            barNetVolumePerMeter = -1;
            barHeight_seaward = -1;
            barWidth_seaward = -1;
            barSlope_seaward = -1;
            barHeight_landward = -1;
            barWidth_landward = -1;
            barSlope_landward = -1;
            ###################################################################
            # CALCULATIONS LANDWARD PART OF INTERTIDAL BAR
            ###################################################################
            # Trough point on the landward part of the crest point has the same ID
            trough_ixCrossSection_landward = self.np_trough_ixCrossSection[
                pointCount
            ];
            trough_topo_landward = self.crossSectionProfile_topo_raw[
                trough_ixCrossSection_landward
            ];
            barHeight_landward = crest_topo - trough_topo_landward;
            barWidth_landward = (trough_ixCrossSection_landward - crest_ixCrossSection)\
                                * self.rPixelResolution_mean;
            barSlope_landward = np.degrees(np.arctan(
                barHeight_landward / barWidth_landward
            ));
            ###################################################################
            # CALCULATIONS SEAWARD PART OF THE INTERTIDAL BAR
            ###################################################################
            # This means the point is the extrema point closest to the sea
            #  Not possible to calculate seaward part of the bar
            if pointCount != 0:
                trough_ID_seaward = self.np_crestTrough_pointID[
                    pointCount - 1
                ];
                trough_ixCrossSection_seaward = self.np_trough_ixCrossSection[
                    pointCount - 1
                ];
                trough_topo_seaward = self.crossSectionProfile_topo_raw[
                    trough_ixCrossSection_seaward
                ];
                
                barHeight_seaward = crest_topo - trough_topo_seaward;
                barWidth_seaward = (crest_ixCrossSection - trough_ixCrossSection_seaward)\
                                   * self.rPixelResolution_mean;
                barSlope_seaward = np.degrees(np.arctan(
                    barHeight_seaward / barWidth_seaward
                ));
                
                barWidth_total = (trough_ixCrossSection_landward - trough_ixCrossSection_seaward)\
                                 * self.rPixelResolution_mean;
                barSymmetry = barWidth_landward / barWidth_total;
                # Net volume per meter is calculated by subtracting the area
                #  under the intertidal bar with the area under the trapezoid
                #  found by using the topography of the trough points before
                #  and after the intertidal bar
                netVolumePerMeter_base = barWidth_total\
                                         * (trough_topo_landward + trough_topo_seaward)/2;
                netVolumePerMeter_bar = np.sum((self.crossSectionProfile_topo_raw[
                    trough_ixCrossSection_seaward:trough_ixCrossSection_landward
                ]) * self.rPixelResolution_mean);
                barNetVolumePerMeter = netVolumePerMeter_bar - netVolumePerMeter_base;
            ###################################################################
            # STORING INFORMATION IN ARRAYS
            ###################################################################
            crossSectionCrest_barCharacteristics = np.array([
                self.raster_crossSection_count, crest_ID,
                np.int(barWidth_total),
                np.around(barSymmetry,3), np.around(barNetVolumePerMeter,2),
                trough_ID_seaward, np.around(barHeight_seaward,2),
                np.int(barWidth_seaward), np.around(barSlope_seaward,3),
                np.around(barHeight_landward,2), np.int(barWidth_landward),
                np.around(barSlope_landward,3)
            ]);
            
            if self.np_raster_barCharacteristics.size == 0:
                self.np_raster_barCharacteristics = crossSectionCrest_barCharacteristics;
            else:
                self.np_raster_barCharacteristics = np.vstack((
                    self.np_raster_barCharacteristics,
                    crossSectionCrest_barCharacteristics
                ));

    def run_through_crossSections_over_raster(
            self,
            pd_lowReferencePoints_on_complete_raster,
            threshold_min_topo_valid,
            min_distance_crest_trough_m,
            np_coordinatesCrossSection_refRaster_start_x,
            np_coordinatesCrossSection_refRaster_start_y,
            np_coordinatesCrossSection_refRaster_end_x,
            np_coordinatesCrossSection_refRaster_end_y,
            freq_plotting_crossSections,
            plot_folderPath_out,
            raster_shape,
            topo_interval):
        """Loop over the available cross-sections for the raster"""
        self.np_raster_crossSection_count = np.arange(
            0, np_coordinatesCrossSection_refRaster_start_x.size, 1,
            dtype=np.int
        );
        # Run through all the available cross-sections
        for self.raster_crossSection_count in self.np_raster_crossSection_count:
            # Find the location (index) of the low reference point on the selected cross-section
            try:
                pd_index_on_crossSectionProfile_ThresholdLower = pd_lowReferencePoints_on_complete_raster.loc[
                    pd_lowReferencePoints_on_complete_raster['CrossSection_ID'] == self.raster_crossSection_count,
                    'Point_CrossSectionIndex'
                ];
            except IndexError:
                index_on_crossSectionProfile_ThresholdLower = 0;
            except KeyError:
                index_on_crossSectionProfile_ThresholdLower = 0;
            else:
                if pd_index_on_crossSectionProfile_ThresholdLower.empty:
                    index_on_crossSectionProfile_ThresholdLower = -1;
                else:
                    index_on_crossSectionProfile_ThresholdLower = np.int(
                        pd_index_on_crossSectionProfile_ThresholdLower.iloc[0]
                    );
            # Extract coordinates of selected cross-section from all available cross-sections.
            self.coordinatesCrossSection_refRaster_start_x = np_coordinatesCrossSection_refRaster_start_x[
                self.raster_crossSection_count
            ];
            self.coordinatesCrossSection_refRaster_start_y = np_coordinatesCrossSection_refRaster_start_y[
                self.raster_crossSection_count
            ];
            self.coordinatesCrossSection_refRaster_end_x = np_coordinatesCrossSection_refRaster_end_x[
                self.raster_crossSection_count
            ];
            self.coordinatesCrossSection_refRaster_end_y = np_coordinatesCrossSection_refRaster_end_y[
                self.raster_crossSection_count
            ];
            
            self.do_analysis_of_crossSectionProfile(
                index_on_crossSectionProfile_ThresholdLower,
                threshold_min_topo_valid,
                min_distance_crest_trough_m,
                self.coordinatesCrossSection_refRaster_start_x,
                self.coordinatesCrossSection_refRaster_start_y,
                self.coordinatesCrossSection_refRaster_end_x,
                self.coordinatesCrossSection_refRaster_end_y
            );
            self._append_inflection_from_crossSection_to_raster();
            self._append_crestTrough_from_crossSection_to_raster();       
            
            if np.mod(self.raster_crossSection_count, freq_plotting_crossSections) != 0\
                    or freq_plotting_crossSections == -1\
                    or self.np_crest_ixCrossSection.size == 0\
                    or plot_folderPath_out is None:
                continue;
            plt_ib.plot_crossSectionProfile_intertidalBar(
                self.raster_topo_raw,
                raster_shape,
                (self.rPixelResolution_x,self.rPixelResolution_y),
                (self.rEasting_origin,self.rNorthing_origin),
                topo_interval,
                (self.coordinatesCrossSection_refRaster_start_x,self.coordinatesCrossSection_refRaster_start_y),
                (self.coordinatesCrossSection_refRaster_end_x,self.coordinatesCrossSection_refRaster_end_y),
                self.raster_crossSection_count,
                self.crossSectionProfile_topo_raw,
                self.crossSectionProfile_topo_smooth,
                self.np_inflection_pointID,
                self.np_inflection_ixCrossSection,
                self.np_crestTrough_pointID,
                self.np_crest_ixCrossSection,
                self.np_crest_x,
                self.np_crest_y,
                self.np_trough_ixCrossSection,
                self.np_trough_x,
                self.np_trough_y,
                plot_folderPath_out
            );
             
    def do_analysis_of_crossSectionProfile(
            self,
            index_on_crossSectionProfile_ThresholdLower,
            threshold_min_topo_valid,
            min_distance_crest_trough_m,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y):
        """Does the analysis for one specific cross-section profile"""
        # Cross-section profile over the raw (non-smoothed) raster
        self.crossSectionProfile_topo_raw,\
        self.rr,\
        self.cc = cs_utils.get_crossSectionProfile_from_raster(
            self.raster_topo_raw,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y
        );
        # Cross-section profile over the smoothed raster
        self.crossSectionProfile_topo_smooth,\
        _, _ = cs_utils.get_crossSectionProfile_from_raster(
            self.raster_topo_smooth,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y
        );
        # Initialises all arrays that contain information from feature points
        self._init_featurePoints_on_crossSectionProfile();
        # Chech whether the length of the cross-section is enough to draw conclusions
        self._isCrossSectionLengthEnough(
            self.crossSectionProfile_topo_smooth
        );
        if self.isCrossSectionLengthEnough:
            # Search the inflection points on the cross-section profile
            self._find_inflectionPoints_on_crossSectionProfile();
            # Searches the crest, trough points on the cross-section profile
            self._find_crestOrTroughPoints_on_crossSectionProfile(
                index_on_crossSectionProfile_ThresholdLower,
                threshold_min_topo_valid,
                min_distance_crest_trough_m
            );
            self._remove_duplicates_localExtremaPoints(True);
            self._barCharacteristics_on_crossSectionProfile();
    
    def get_crossSections_over_complete_raster(
            self,
            np_coordinatesCrossSection_refRaster_start_x,
            np_coordinatesCrossSection_refRaster_start_y,
            np_coordinatesCrossSection_refRaster_end_x,
            np_coordinatesCrossSection_refRaster_end_y):
        
        np_coordinatesCrossSection_refEN_start_easting = self.rEasting_origin\
            + np_coordinatesCrossSection_refRaster_start_x*self.rPixelResolution_x;
        np_coordinatesCrossSection_refEN_start_northing = self.rNorthing_origin\
            + np_coordinatesCrossSection_refRaster_start_y*self.rPixelResolution_y;
        np_coordinatesCrossSection_refEN_end_easting = self.rEasting_origin\
            + np_coordinatesCrossSection_refRaster_end_x*self.rPixelResolution_x;
        np_coordinatesCrossSection_refEN_end_northing = self.rNorthing_origin\
            + np_coordinatesCrossSection_refRaster_end_y*self.rPixelResolution_y;
        crossSections_over_complete_raster = np.column_stack((
            self.np_raster_crossSection_count,
            np_coordinatesCrossSection_refRaster_start_x,
            np_coordinatesCrossSection_refRaster_start_y,
            np_coordinatesCrossSection_refEN_start_easting,
            np_coordinatesCrossSection_refEN_start_northing,
            np_coordinatesCrossSection_refRaster_end_x,
            np_coordinatesCrossSection_refRaster_end_y,
            np_coordinatesCrossSection_refEN_end_easting,
            np_coordinatesCrossSection_refEN_end_northing
        ));
        colsHeader = np.array([
            'CrossSection_ID', 'CrossSection_StartX', 'CrossSection_StartY',
            'CrossSection_StartEasting', 'CrossSection_StartNorthing',
            'CrossSection_EndX', 'CrossSection_EndY',
            'CrossSection_EndEasting', 'CrossSection_EndNorthing'
        ]);
        pd_return = pd.DataFrame(
            data=crossSections_over_complete_raster,
            columns=colsHeader
        );
        return pd_return;
    
    def get_inflectionPoints_over_complete_raster(self):
        np_raster_inflection_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_inflection_x;
        np_raster_inflection_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_inflection_y;
        inflectionPoints_over_complete_raster = np.column_stack((
            self.np_raster_inflection_crossSection_ID,
            self.np_raster_inflection_ID, self.np_raster_inflection_ixCrossSection,
            self.np_raster_inflection_x, self.np_raster_inflection_y,
            np_raster_inflection_easting, np_raster_inflection_northing,
            self.np_raster_inflection_topo
        ));
        colsHeader = np.array([
            'CrossSection_ID',
            'Point_ID', 'Point_CrossSectionIndex',
            'Point_X', 'Point_Y',
            'Point_Easting', 'Point_Northing',
            'Point_Topography'
        ]);
        pd_return = pd.DataFrame(
            data=inflectionPoints_over_complete_raster,
            columns=colsHeader
        );
        return pd_return;
    
    def get_crestPoints_over_complete_raster(self):
        np_raster_crest_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_crest_x;
        np_raster_crest_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_crest_y;        
        crestPoints_over_complete_raster = np.column_stack((
            self.np_raster_crestTrough_crossSection_ID,
            self.np_raster_crestTrough_ID, self.np_raster_crest_ixCrossSection,
            self.np_raster_crest_x, self.np_raster_crest_y,
            np_raster_crest_easting, np_raster_crest_northing,
            self.np_raster_crest_topo
        ));
        colsHeader = np.array([
            'CrossSection_ID',
            'Point_ID','Point_CrossSectionIndex',
            'Point_X','Point_Y',
            'Point_Easting','Point_Northing',
            'Point_Topography'
        ]);
        return pd.DataFrame(
            data=crestPoints_over_complete_raster,
            columns=colsHeader
        );
        
    def get_troughPoints_over_complete_raster(self):
        np_raster_trough_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_trough_x;
        np_raster_trough_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_trough_y;        
        troughPoints_over_complete_raster = np.column_stack((
            self.np_raster_crestTrough_crossSection_ID,
            self.np_raster_crestTrough_ID, self.np_raster_trough_ixCrossSection,
            self.np_raster_trough_x, self.np_raster_trough_y,
            np_raster_trough_easting, np_raster_trough_northing,
            self.np_raster_trough_topo
        ));
        colsHeader = np.array([
            'CrossSection_ID',
            'Point_ID','Point_CrossSectionIndex',
            'Point_X','Point_Y',
            'Point_Easting','Point_Northing',
            'Point_Topography'
        ]);

        return pd.DataFrame(
            data=troughPoints_over_complete_raster,
            columns=colsHeader
        );

    def get_barCharacteristics_over_complete_raster(self):
        colHeaders = np.array([
            'CrossSection_ID', 'Crest_PointID',
            'Bar_Width_Total', 'Bar_Symmetry', 'Bar_CrossSection_NetVolumePerMeter',
            'Trough_PointID_Seaward', 'Bar_Height_Seaward', 'Bar_Width_Seaward',
            'Bar_Slope_Seaward', 'Bar_Height_Landward', 'Bar_Width_Landward',
            'Bar_Slope_Landward'
        ]);
        return pd.DataFrame(
            data=self.np_raster_barCharacteristics,
            columns=colHeaders
        );
