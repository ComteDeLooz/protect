# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy import signal

import config
import cross_section_analysis.utils as cs_utils
import constants
import plotting.plot_drainage_channel as plt_dc


class ChannelPointsOnRaster():
    """
    Taking cross-sections parallel to the beach in order to identify
     drainage channels cutting through the intertidal bars
    """
    def __init__(
            self,
            raster_pixel_resolution, raster_origin_refEN,
            raster_topo_raw, raster_topo_smooth,
            minimum_channel_width_m, minimum_channel_depth_m
            ):
        
        self.rPixelResolution_x,\
        self.rPixelResolution_y = raster_pixel_resolution;
        self.rPixelResolution_mean = np.mean(
            np.absolute([self.rPixelResolution_x, self.rPixelResolution_y])
        );
        self.rEasting_origin, self.rNorthing_origin = raster_origin_refEN;
        self.raster_topo_raw = raster_topo_raw;
        self.raster_topo_smooth = raster_topo_smooth;
        self.minimum_channel_width_pxl = minimum_channel_width_m\
                                         / self.rPixelResolution_mean;
        self.minimum_channel_depth_m = minimum_channel_depth_m;
        self.init_featurePoints_on_raster();

    def _isCrossSectionLengthEnough(self, crossSectionProfile_topo):
        """If a cross-section is too short, no relevant data are obtained"""
        minLength_pxl = config.CROSS_SECTION_SEARCH_CHANNEL_FEATURES['MINIMUM_LENGTH_CROSSSECTION_METER']\
                        / self.rPixelResolution_mean;
        self.isCrossSectionLengthEnough = True;
        if crossSectionProfile_topo.size <= minLength_pxl:
            self.isCrossSectionLengthEnough = False;

    def _init_featurePoints_on_crossSectionProfile(self):
        """
        Initialisation of new arrays that contain channel features found over the
         selected cross-sections
        Re-initializes for each new cross-section
        """
        # All inflection points found, also those not linked to a drainage channel
        self.np_inflectionGeneral_pointID = np.array([], dtype=np.int);
        self.np_inflectionGeneral_x = np.array([], dtype=np.int);
        self.np_inflectionGeneral_y = np.array([], dtype=np.int);
        self.np_inflectionLeft_topo = np.array([], dtype=np.float);
        self.np_inflectionGeneral_ixCrossSection = np.array([], dtype=np.int);

        self.np_channel_pointID = np.array([], dtype=np.int);
        
        # Inflection points thought to be at the left slope of the channels
        #  (looking towards the sea)
        self.np_inflectionLeft_x = np.array([], dtype=np.int);
        self.np_inflectionLeft_y = np.array([], dtype=np.int);
        self.np_inflectionLeft_topo = np.array([], dtype=np.float);
        self.np_inflectionLeft_ixCrossSection = np.array([], dtype=np.int);
        
        # Inflection points thought to be at the right slope of the channels
        #  (looking towards the sea)
        self.np_inflectionRight_x = np.array([], dtype=np.int);
        self.np_inflectionRight_y = np.array([], dtype=np.int);
        self.np_inflectionRight_topo = np.array([], dtype=np.float);
        self.np_inflectionRight_ixCrossSection = np.array([], dtype=np.int);
        
        # Points at the bottom of the channels
        self.np_bottom_x = np.array([], dtype=np.int);
        self.np_bottom_y = np.array([], dtype=np.int);
        self.np_bottom_topo = np.array([], dtype=np.float);
        self.np_bottom_ixCrossSection = np.array([], dtype=np.int);
        
        # Points at the top left of the slope of the channels (left end)
        #  (looking towards the sea)
        self.np_topLeft_x = np.array([], dtype=np.int);
        self.np_topLeft_y = np.array([], dtype=np.int);
        self.np_topLeft_topo = np.array([], dtype=np.float);
        self.np_topLeft_ixCrossSection = np.array([], dtype=np.int);
        
        # Points at the top right of the slope of the channels (right end)
        #  (looking towards the sea)
        self.np_topRight_x = np.array([], dtype=np.int);
        self.np_topRight_y = np.array([], dtype=np.int);
        self.np_topRight_topo = np.array([], dtype=np.float);
        self.np_topRight_ixCrossSection = np.array([], dtype=np.int);

    def init_featurePoints_on_raster(self):
        """
        Initialisation of numpy arrays that store all channels feature information
         over the complete topography raster
        """
        self.featurePoint_ID = 900000;
        
        self.np_raster_inflectionGeneral_pointID = np.array([], dtype=np.int);
        self.np_raster_inflectionGeneral_x = np.array([], dtype=np.int);
        self.np_raster_inflectionGeneral_y = np.array([], dtype=np.int);
        self.np_raster_inflectionGeneral_topo = np.array([], dtype=np.float);
        self.np_raster_inflectionGeneral_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_inflectionGeneral_crossSection_ID = np.array([], dtype=np.int);
        
        self.np_raster_channel_pointID = np.array([], dtype=np.int);
        self.np_raster_channel_crossSection_ID = np.array([], dtype=np.int);
        
        self.np_raster_inflectionLeft_x = np.array([], dtype=np.int);
        self.np_raster_inflectionLeft_y = np.array([], dtype=np.int);
        self.np_raster_inflectionLeft_topo = np.array([], dtype=np.float);
        self.np_raster_inflectionLeft_ixCrossSection = np.array([], dtype=np.int);
        
        self.np_raster_inflectionRight_x = np.array([], dtype=np.int);
        self.np_raster_inflectionRight_y = np.array([], dtype=np.int);
        self.np_raster_inflectionRight_topo = np.array([], dtype=np.float);
        self.np_raster_inflectionRight_ixCrossSection = np.array([], dtype=np.int);
        
        self.np_raster_bottom_raw_x = np.array([], dtype=np.int);
        self.np_raster_bottom_raw_y = np.array([], dtype=np.int);
        self.np_raster_bottom_raw_topo = np.array([], dtype=np.float);
        self.np_raster_bottom_raw_ixCrossSection = np.array([], dtype=np.int);
        
        self.np_raster_topLeft_x = np.array([], dtype=np.int);
        self.np_raster_topLeft_y = np.array([], dtype=np.int);
        self.np_raster_topLeft_topo = np.array([], dtype=np.float);
        self.np_raster_topLeft_ixCrossSection = np.array([], dtype=np.int);
        
        self.np_raster_topRight_x = np.array([], dtype=np.int);
        self.np_raster_topRight_y = np.array([], dtype=np.int);
        self.np_raster_topRight_topo = np.array([], dtype=np.float);
        self.np_raster_topRight_ixCrossSection = np.array([], dtype=np.int);

    def _find_inflectionPoints_on_crossSectionProfile(self):
        """Searches inflection points on the cross-section profile."""        
        profile_topo_firstDeriv = np.gradient(
            self.crossSectionProfile_topo_smooth
        );
        profile_topo_secondDeriv = np.gradient(
            profile_topo_firstDeriv
        );
        
        medianFilter_kernelSize_pxl = np.int(np.around(
            constants.CHANNEL['INFLECTION_POINTS_MEDIAN_FILTER_KERNEL_METER']\
            / self.rPixelResolution_mean)
        );
        if np.mod(medianFilter_kernelSize_pxl, 2) == 0:
            medianFilter_kernelSize_pxl += 1;
        profile_topo_secondDeriv_smooth = signal.medfilt(
            profile_topo_secondDeriv, kernel_size=medianFilter_kernelSize_pxl
        );
        inflectionGeneral_ixCrossSection = np.where(
            np.diff(np.sign(profile_topo_secondDeriv_smooth)) != 0
        )[0];
        # Stores all indices of inflection points found on the cross-section
        # but only for the considered cross-section. If another cross-section
        # is loaded, this array re-initializes.
        if inflectionGeneral_ixCrossSection.size != 0:
            self.np_inflectionGeneral_ixCrossSection = np.sort(
                inflectionGeneral_ixCrossSection
            );     
    
    def _append_inflection_from_crossSection_to_raster(self):
        """
        Adds info on the inflection points found on the considered
        cross-section to the array containing all inflection points for the
        complete raster.
        """
        self.np_raster_inflectionGeneral_pointID = np.append(
            self.np_raster_inflectionGeneral_pointID,
            self.np_inflectionGeneral_pointID
        );
        self.np_raster_inflectionGeneral_x = np.append(
            self.np_raster_inflectionGeneral_x,
            self.cc[self.np_inflectionGeneral_ixCrossSection]
        );
        self.np_raster_inflectionGeneral_y = np.append(
            self.np_raster_inflectionGeneral_y,
            self.rr[self.np_inflectionGeneral_ixCrossSection]
        );
        self.np_raster_inflectionGeneral_topo = np.append(
            self.np_raster_inflectionGeneral_topo,
            self.crossSectionProfile_topo_raw[self.np_inflectionGeneral_ixCrossSection]
        );
        self.np_raster_inflectionGeneral_ixCrossSection = np.append(
            self.np_raster_inflectionGeneral_ixCrossSection,
            self.np_inflectionGeneral_ixCrossSection
        );
        
        numberOfInflectionPoints = self.np_inflectionGeneral_ixCrossSection.size;
        self.np_raster_inflectionGeneral_crossSection_ID = np.append(
            self.np_raster_inflectionGeneral_crossSection_ID,
            np.full((numberOfInflectionPoints,), self.raster_crossSection_count, dtype=np.int)
        );

    def _detect_channels_from_inflection_points_on_crossSectionProfile(self):
        """
        A channel is considered to have the following features:
            * Top of the channel (2x)
                - tl = top left of the channel
                - tr = top right of the channel
            * Bottom of the channel (lowest topo) (1x)
                - b = bottom of the channel
            * Inflection on the slope between top and bottom (2x)
                - il = left inflection of the channel
                - ir = right inflection of the channel
        
        To find the features above, the inflection point array is used,
        deduced from the method _find_all_inflection_points_on_crossSection().
        
        The following inflection points are used to detect a channel from a 
        topography profile (see figure below):
            * i1 = start channel search
            * i2 = left inflection point of the channel
            * i3 = right inflection point of the channel
            * i4 = end channel search
        
                                  --tr--
                                 /      \
                 --tl---        /        i4    /
                /       \      ir=i3      \   /
               /     i2=il    /            ---
         \    i1         \   /
          \  /            -b-
           --
        """
        # Loop through the inflection points found for the considered cross-section.
        for (index,), i2_ixCrossSection in np.ndenumerate(
                self.np_inflectionGeneral_ixCrossSection):
            # Unique ID identifier per inflection point found on the raster.
            crossSection_inflection_pointID = self.featurePoint_ID;
            self.np_inflectionGeneral_pointID = np.append(
                self.np_inflectionGeneral_pointID,
                crossSection_inflection_pointID
            );
            self.featurePoint_ID += 1;
            # The first inflection point on the cross-section is not considered
            # to be the left side of a channel.
            if index == 0:
                continue;
            try:
                # Coordinate on the cross-section of the second next inflection point
                i4_ixCrossSection = self.np_inflectionGeneral_ixCrossSection[
                    index + 2
                ];
            # When too close to the end of the cross-section.
            except IndexError:
                continue;       
            # Coordinate on the cross-section of the previous inflection point
            i1_ixCrossSection = self.np_inflectionGeneral_ixCrossSection[
                index - 1
            ];
            # Coordinate on the cross-section of the next inflection point
            i3_ixCrossSection = self.np_inflectionGeneral_ixCrossSection[
                index + 1
            ];
            
            try:
                # position of bottom of channel based on non-smoothed raster
                b_ixCrossSection = i2_ixCrossSection + np.argmin(
                    self.crossSectionProfile_topo_raw[i2_ixCrossSection:i3_ixCrossSection]
                );
                # position of top left of channel on non-smoothed raster
                tl_ixCrossSection = i1_ixCrossSection + np.argmax(
                    self.crossSectionProfile_topo_raw[i1_ixCrossSection:i2_ixCrossSection]
                );
                # position of top right of channel on non-smoothed raster
                tr_ixCrossSection = i3_ixCrossSection + np.argmax(
                    self.crossSectionProfile_topo_raw[i3_ixCrossSection:i4_ixCrossSection]
                );
            except ValueError:
                continue;
            # Condition 1:
            # Position of bottom, top left and top right of the channel should
            #  be different from position of inflection points
            if (b_ixCrossSection == i2_ixCrossSection)\
                    or (b_ixCrossSection == i3_ixCrossSection)\
                    or (tl_ixCrossSection == i2_ixCrossSection)\
                    or (tr_ixCrossSection == i3_ixCrossSection):
                continue;
            
            i2_topo = self.crossSectionProfile_topo_raw[i2_ixCrossSection];
            i3_topo = self.crossSectionProfile_topo_raw[i3_ixCrossSection];
            b_topo = self.crossSectionProfile_topo_raw[b_ixCrossSection];
            tl_topo = self.crossSectionProfile_topo_raw[tl_ixCrossSection];
            tr_topo = self.crossSectionProfile_topo_raw[tr_ixCrossSection];
            # Condition 2: 
            # Topography of inflection point on the left side should be lower
            #  than the topography of the maximum of this same channel side.
            # Similar for the right side of the channel.
            if (i2_topo >= tl_topo) or (i3_topo >= tr_topo):
                continue;
            # Condition 3:
            # Difference between top and bottom of channel should be at least
            #  the specified height
            if (tl_topo-b_topo) <= self.minimum_channel_depth_m\
                    or (tr_topo-b_topo) <= self.minimum_channel_depth_m:
                continue;
            # Condition 4:
            # Width of the channel (measured at the top) should be at least
            #  the user-specified width
            if (tr_ixCrossSection-tl_ixCrossSection) <= self.minimum_channel_width_pxl:
                continue;
            
            # Store information in arrays related to the cross-section
            self.np_channel_pointID = np.append(
                self.np_channel_pointID, crossSection_inflection_pointID
            );
            self.np_inflectionLeft_ixCrossSection = np.append(
                self.np_inflectionLeft_ixCrossSection, i2_ixCrossSection
            );
            self.np_inflectionRight_ixCrossSection = np.append(
                self.np_inflectionRight_ixCrossSection, i3_ixCrossSection
            );
            self.np_bottom_ixCrossSection = np.append(
                self.np_bottom_ixCrossSection, b_ixCrossSection
            );
            self.np_topLeft_ixCrossSection = np.append(
                self.np_topLeft_ixCrossSection, tl_ixCrossSection
            );
            self.np_topRight_ixCrossSection = np.append(
                self.np_topRight_ixCrossSection, tr_ixCrossSection
            );

    def _append_channel_from_crossSection_to_raster(self):
        """
        Adds info on the channels found on the considered cross-section
        to the array containing all inflection points for the complete raster.
        """        
        self.np_raster_channel_pointID = np.append(
            self.np_raster_channel_pointID,
            self.np_channel_pointID);
        
        self.np_raster_inflectionLeft_x = np.append(
            self.np_raster_inflectionLeft_x,
            self.cc[self.np_inflectionLeft_ixCrossSection]
        );
        self.np_raster_inflectionLeft_y = np.append(
            self.np_raster_inflectionLeft_y,
            self.rr[self.np_inflectionLeft_ixCrossSection]
        );
        self.np_raster_inflectionLeft_topo = np.append(
            self.np_raster_inflectionLeft_topo,
            self.crossSectionProfile_topo_raw[self.np_inflectionLeft_ixCrossSection]
        );
        self.np_raster_inflectionLeft_ixCrossSection = np.append(
            self.np_raster_inflectionLeft_ixCrossSection,
            self.np_inflectionLeft_ixCrossSection
        );
        
        self.np_raster_inflectionRight_x = np.append(
            self.np_raster_inflectionRight_x,
            self.cc[self.np_inflectionRight_ixCrossSection]
        );
        self.np_raster_inflectionRight_y = np.append(
            self.np_raster_inflectionRight_y,
            self.rr[self.np_inflectionRight_ixCrossSection]
        );
        self.np_raster_inflectionRight_topo = np.append(
            self.np_raster_inflectionRight_topo,
            self.crossSectionProfile_topo_raw[self.np_inflectionRight_ixCrossSection]
        );
        self.np_raster_inflectionRight_ixCrossSection = np.append(
            self.np_raster_inflectionRight_ixCrossSection,
            self.np_inflectionRight_ixCrossSection
        );
        
        self.np_raster_bottom_raw_x = np.append(
            self.np_raster_bottom_raw_x,
            self.cc[self.np_bottom_ixCrossSection]
        );
        self.np_raster_bottom_raw_y = np.append(
            self.np_raster_bottom_raw_y,
            self.rr[self.np_bottom_ixCrossSection]
        );
        self.np_raster_bottom_raw_topo = np.append(
            self.np_raster_bottom_raw_topo,
            self.crossSectionProfile_topo_raw[self.np_bottom_ixCrossSection]
        );
        self.np_raster_bottom_raw_ixCrossSection = np.append(
            self.np_raster_bottom_raw_ixCrossSection,
            self.np_bottom_ixCrossSection
        );
                
        self.np_raster_topLeft_x = np.append(
            self.np_raster_topLeft_x,
            self.cc[self.np_topLeft_ixCrossSection]
        );
        self.np_raster_topLeft_y = np.append(
            self.np_raster_topLeft_y,
            self.rr[self.np_topLeft_ixCrossSection]
        );
        self.np_raster_topLeft_topo = np.append(
            self.np_raster_topLeft_topo,
            self.crossSectionProfile_topo_raw[self.np_topLeft_ixCrossSection]
        );
        self.np_raster_topLeft_ixCrossSection = np.append(
            self.np_raster_topLeft_ixCrossSection,
            self.np_topLeft_ixCrossSection
        );
        
        self.np_raster_topRight_x = np.append(
            self.np_raster_topRight_x,
            self.cc[self.np_topRight_ixCrossSection]
        );
        self.np_raster_topRight_y = np.append(
            self.np_raster_topRight_y,
            self.rr[self.np_topRight_ixCrossSection]
        );
        self.np_raster_topRight_topo = np.append(
            self.np_raster_topRight_topo,
            self.crossSectionProfile_topo_raw[self.np_topRight_ixCrossSection]
        );
        self.np_raster_topRight_ixCrossSection = np.append(
            self.np_raster_topRight_ixCrossSection,
            self.np_topRight_ixCrossSection
        );
        
        numberOfChannelPoints = self.np_inflectionLeft_ixCrossSection.size;
        self.np_raster_channel_crossSection_ID = np.append(
            self.np_raster_channel_crossSection_ID,
            np.full((numberOfChannelPoints,), self.raster_crossSection_count, dtype=np.int)
        );

    def run_through_crossSections_over_raster(
            self,
            np_coordinatesCrossSection_refRaster_start_x,
            np_coordinatesCrossSection_refRaster_start_y,
            np_coordinatesCrossSection_refRaster_end_x,
            np_coordinatesCrossSection_refRaster_end_y,
            freq_plotting_crossSections,
            plot_folderPath_out,
            raster_shape,
            topo_interval):
        self.np_raster_crossSection_count = np.arange(
            0, np_coordinatesCrossSection_refRaster_start_x.size, 1,
            dtype=np.int
        );
        for self.raster_crossSection_count in self.np_raster_crossSection_count:
            # Extract coordinates of selected cross-section from all available
            # cross-sections.
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
                self.coordinatesCrossSection_refRaster_start_x,
                self.coordinatesCrossSection_refRaster_start_y,
                self.coordinatesCrossSection_refRaster_end_x,
                self.coordinatesCrossSection_refRaster_end_y
            );
            self._append_inflection_from_crossSection_to_raster();
            self._append_channel_from_crossSection_to_raster();
            
            if np.mod(self.raster_crossSection_count, freq_plotting_crossSections) != 0\
                    or freq_plotting_crossSections == -1\
                    or self.np_bottom_ixCrossSection.size == 0\
                    or plot_folderPath_out is None:
                continue;
            plt_dc.plot_crossSectionProfile_drainageChannel(
                self.raster_topo_raw,
                raster_shape,
                (self.rPixelResolution_x,self.rPixelResolution_y),
                (self.rEasting_origin,self.rNorthing_origin),
                topo_interval[0],
                topo_interval[1],
                (self.coordinatesCrossSection_refRaster_start_x,self.coordinatesCrossSection_refRaster_start_y),
                (self.coordinatesCrossSection_refRaster_end_x,self.coordinatesCrossSection_refRaster_end_y),
                self.raster_crossSection_count,
                self.crossSectionProfile_topo_raw,
                self.crossSectionProfile_topo_smooth,
                self.np_inflectionGeneral_ixCrossSection,
                self.np_bottom_ixCrossSection,
                self.np_raster_bottom_raw_x,
                self.np_raster_bottom_raw_y,
                self.np_topLeft_ixCrossSection,
                self.np_raster_topLeft_x,
                self.np_raster_topLeft_y,
                self.np_topRight_ixCrossSection,
                self.np_raster_topRight_x,
                self.np_raster_topRight_y,
                plot_folderPath_out
            );
    def do_analysis_of_crossSectionProfile(
            self,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y):
        
        self.crossSectionProfile_topo_raw,\
        self.rr,\
        self.cc = cs_utils.get_crossSectionProfile_from_raster(
            self.raster_topo_raw,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y
        );
        self.crossSectionProfile_topo_smooth,\
        _, _ = cs_utils.get_crossSectionProfile_from_raster(
            self.raster_topo_smooth,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y
        );
        self._init_featurePoints_on_crossSectionProfile();
        self._isCrossSectionLengthEnough(
            self.crossSectionProfile_topo_smooth
        );
        if self.isCrossSectionLengthEnough:
            self._find_inflectionPoints_on_crossSectionProfile();
            self._detect_channels_from_inflection_points_on_crossSectionProfile();

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
        """Outputs ass info related to the inflection points for the complete raster.""" 
        np_raster_inflectionGeneral_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_inflectionGeneral_x;
        np_raster_inflectionGeneral_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_inflectionGeneral_y;
        inflectionGeneral = np.column_stack((
            self.np_raster_inflectionGeneral_crossSection_ID,
            self.np_raster_inflectionGeneral_pointID,
            self.np_raster_inflectionGeneral_ixCrossSection,
            self.np_raster_inflectionGeneral_x,
            self.np_raster_inflectionGeneral_y,
            np_raster_inflectionGeneral_easting,
            np_raster_inflectionGeneral_northing,
            self.np_raster_inflectionGeneral_topo
        ));
        colsHeader = np.array([
            'CrossSection_ID',
            'Point_ID','Point_CrossSectionIndex',
            'Point_X', 'Point_Y',
            'Point_Easting', 'Point_Northing',
            'Point_Topography'
        ]);
        pd_return = pd.DataFrame(
            data=inflectionGeneral,
            columns=colsHeader
        );
        return pd_return;

    def get_channelPoints_over_complete_raster(self):
        """Outputs all info related to the channels for the complete raster."""
        np_raster_inflectionLeft_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_inflectionLeft_x;
        np_raster_inflectionLeft_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_inflectionLeft_y;
        inflectionLeftPoints_over_complete_raster = np.column_stack((
            self.np_raster_channel_crossSection_ID,
            self.np_raster_channel_pointID,
            self.np_raster_inflectionLeft_ixCrossSection,
            self.np_raster_inflectionLeft_x, self.np_raster_inflectionLeft_y,
            np_raster_inflectionLeft_easting, np_raster_inflectionLeft_northing,
            self.np_raster_inflectionLeft_topo
        ));
        
        np_raster_inflectionRight_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_inflectionRight_x;
        np_raster_inflectionRight_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_inflectionRight_y;
        inflectionRightPoints_over_complete_raster = np.column_stack((
            self.np_raster_channel_crossSection_ID,
            self.np_raster_channel_pointID,
            self.np_raster_inflectionRight_ixCrossSection,
            self.np_raster_inflectionRight_x, self.np_raster_inflectionRight_y,
            np_raster_inflectionRight_easting, np_raster_inflectionRight_northing,
            self.np_raster_inflectionRight_topo
        ));                                  
        
        np_raster_bottom_raw_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_bottom_raw_x;
        np_raster_bottom_raw_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_bottom_raw_y;
        bottomPoints_raw_over_complete_raster = np.column_stack((
            self.np_raster_channel_crossSection_ID,
            self.np_raster_channel_pointID,
            self.np_raster_bottom_raw_ixCrossSection,
            self.np_raster_bottom_raw_x, self.np_raster_bottom_raw_y,
            np_raster_bottom_raw_easting, np_raster_bottom_raw_northing,
            self.np_raster_bottom_raw_topo
        ));
                
        np_raster_topLeft_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_topLeft_x;
        np_raster_topLeft_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_topLeft_y;
        topLeftPoints_over_complete_raster = np.column_stack((
            self.np_raster_channel_crossSection_ID,
            self.np_raster_channel_pointID,
            self.np_raster_topLeft_ixCrossSection,
            self.np_raster_topLeft_x, self.np_raster_topLeft_y,
            np_raster_topLeft_easting, np_raster_topLeft_northing,
            self.np_raster_topLeft_topo
        ));
        
        np_raster_topRight_easting = self.rEasting_origin\
            + self.rPixelResolution_x*self.np_raster_topRight_x;
        np_raster_topRight_northing = self.rNorthing_origin\
            + self.rPixelResolution_y*self.np_raster_topRight_y;
        topRightPoints_over_complete_raster = np.column_stack((
            self.np_raster_channel_crossSection_ID,
            self.np_raster_channel_pointID,
            self.np_raster_topRight_ixCrossSection,
            self.np_raster_topRight_x, self.np_raster_topRight_y,
            np_raster_topRight_easting, np_raster_topRight_northing,
            self.np_raster_topRight_topo
        ));
        
        colsHeader = np.array([
            'CrossSection_ID',
            'Point_ID',
            'Point_CrossSectionIndex',
            'Point_X', 'Point_Y',
            'Point_Easting', 'Point_Northing',
            'Point_Topography'
        ]);
        pd_il = pd.DataFrame(
            data=inflectionLeftPoints_over_complete_raster,
            columns=colsHeader
        );
        pd_ir = pd.DataFrame(
            data=inflectionRightPoints_over_complete_raster,
            columns=colsHeader
        );
        pd_b_raw = pd.DataFrame(
            data=bottomPoints_raw_over_complete_raster,
            columns=colsHeader
        );
        pd_tl = pd.DataFrame(
            data=topLeftPoints_over_complete_raster,
            columns=colsHeader
        );
        pd_tr = pd.DataFrame(
            data=topRightPoints_over_complete_raster,
            columns=colsHeader
        );
        return pd_il, pd_ir, pd_b_raw, pd_tl, pd_tr;
