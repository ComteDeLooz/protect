# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import cross_section_analysis.utils as cs_utils


class WaterReferenceLinesOnRaster():
    """
    Processing of cross-section profiles to find water reference lines on the
     topography raster within user-specified topography intervals:
         * Low reference line
         * mean low water line
         * mean intertidal line
         * mean high water line
         * highest astronomical tide line
    """
    def __init__(
            self,
            raster_pixel_resolution,
            raster_origin_refEN,
            raster_topo_raw,
            searchInterval_lowReferenceLine=(),
            searchInterval_meanLowWater=(),
            searchInterval_meanIntertidal=(),
            searchInterval_meanHighWater=(),
            searchInterval_highestAstronomicalTide=()):
        
        self.rPixelResolution_x, self.rPixelResolution_y = raster_pixel_resolution;
        self.rPixelResolution_mean = np.mean(
            np.absolute([self.rPixelResolution_x, self.rPixelResolution_y])
        );
        self.rEasting_origin, self.rNorthing_origin = raster_origin_refEN;
        self.raster_topo_raw = raster_topo_raw;
        
        self.searchInterval_lowReferenceLine = searchInterval_lowReferenceLine;
        self.searchInterval_meanLowWater = searchInterval_meanLowWater;
        self.searchInterval_meanIntertidal = searchInterval_meanIntertidal;
        self.searchInterval_meanHighWater = searchInterval_meanHighWater;
        self.searchInterval_highestAstronomicalTide = searchInterval_highestAstronomicalTide;

    def _init_markerPoints_on_crossSectionProfile(self):
        # Low Reference Line
        # Location (index) on the cross-section of the low reference points
        self.np_raster_lowReferenceLine_ixCrossSection = np.array([], dtype=np.int);
        # Column numbers (X) of the image containing the low reference points
        self.np_raster_lowReferenceLine_x = np.array([], dtype=np.int);
        # Row numbers (Y) of the image containing the low reference points
        self.np_raster_lowReferenceLine_y = np.array([], dtype=np.int);
        # Topography values at the location of the low reference points
        self.np_raster_lowReferenceLine_topo = np.array([], dtype=np.float32);
        # ID of the cross-section where the low reference points are found on
        self.np_raster_lowReferenceLine_crossSection_ID = np.array([], dtype=np.int);
        # Mean Low Water
        self.np_raster_meanLowWater_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_meanLowWater_x = np.array([], dtype=np.int);
        self.np_raster_meanLowWater_y = np.array([], dtype=np.int);
        self.np_raster_meanLowWater_topo = np.array([], dtype=np.float32);
        self.np_raster_meanLowWater_crossSection_ID = np.array([], dtype=np.int);
        # Mean Intertidal Line
        self.np_raster_meanIntertidal_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_meanIntertidal_x = np.array([], dtype=np.int);
        self.np_raster_meanIntertidal_y = np.array([], dtype=np.int);
        self.np_raster_meanIntertidal_topo = np.array([], dtype=np.float32);
        self.np_raster_meanIntertidal_crossSection_ID = np.array([], dtype=np.int);
        # Mean High Water
        self.np_raster_meanHighWater_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_meanHighWater_x = np.array([], dtype=np.int);
        self.np_raster_meanHighWater_y = np.array([], dtype=np.int);
        self.np_raster_meanHighWater_topo = np.array([], dtype=np.float32);
        self.np_raster_meanHighWater_crossSection_ID = np.array([], dtype=np.int);
        # Highest Astronomical Tide
        self.np_raster_highestAstronomicalTide_ixCrossSection = np.array([], dtype=np.int);
        self.np_raster_highestAstronomicalTide_x = np.array([], dtype=np.int);
        self.np_raster_highestAstronomicalTide_y = np.array([], dtype=np.int);
        self.np_raster_highestAstronomicalTide_topo = np.array([], dtype=np.float32);
        self.np_raster_highestAstronomicalTide_crossSection_ID = np.array([], dtype=np.int);
        
        self.marker_colsHeader = [
                'CrossSection_ID', 'Point_CrossSectionIndex',
                'Point_X', 'Point_Y',
                'Point_Easting', 'Point_Northing',
                'Point_Topography'
            ];

    def _append_lowReferencePoint_from_crossSection_to_raster(
            self, point_ixCrossSection):
        point_x = self.cc[point_ixCrossSection];
        point_y = self.rr[point_ixCrossSection];
        self.np_raster_lowReferenceLine_x = np.append(
            self.np_raster_lowReferenceLine_x,
            point_x
        );
        self.np_raster_lowReferenceLine_y = np.append(
            self.np_raster_lowReferenceLine_y,
            point_y
        );
        self.np_raster_lowReferenceLine_topo = np.append(
            self.np_raster_lowReferenceLine_topo,
            self.raster_topo_raw[point_y, point_x]
        );
        self.np_raster_lowReferenceLine_ixCrossSection = np.append(
            self.np_raster_lowReferenceLine_ixCrossSection,
            point_ixCrossSection
        );
        self.np_raster_lowReferenceLine_crossSection_ID = np.append(
            self.np_raster_lowReferenceLine_crossSection_ID,
            self.raster_crossSection_count
        );
    
    def _append_meanLowWaterPoint_from_crossSection_to_raster(
            self, point_ixCrossSection):
        point_x = self.cc[point_ixCrossSection];
        point_y = self.rr[point_ixCrossSection];
        self.np_raster_meanLowWater_x = np.append(
            self.np_raster_meanLowWater_x,
            point_x
        );
        self.np_raster_meanLowWater_y = np.append(
            self.np_raster_meanLowWater_y,
            point_y
        );
        self.np_raster_meanLowWater_topo = np.append(
            self.np_raster_meanLowWater_topo,
            self.raster_topo_raw[point_y, point_x]
        );
        self.np_raster_meanLowWater_ixCrossSection = np.append(
            self.np_raster_meanLowWater_ixCrossSection,
            point_ixCrossSection
        );
        self.np_raster_meanLowWater_crossSection_ID = np.append(
            self.np_raster_meanLowWater_crossSection_ID,
            self.raster_crossSection_count
        );

    def _append_meanIntertidalPoint_from_crossSection_to_raster(
            self, point_ixCrossSection):
        point_x = self.cc[point_ixCrossSection];
        point_y = self.rr[point_ixCrossSection];
        self.np_raster_meanIntertidal_x = np.append(
            self.np_raster_meanIntertidal_x,
            point_x
        );
        self.np_raster_meanIntertidal_y = np.append(
            self.np_raster_meanIntertidal_y,
            point_y
        );
        self.np_raster_meanIntertidal_topo = np.append(
            self.np_raster_meanIntertidal_topo,
            self.raster_topo_raw[point_y, point_x]
        );
        self.np_raster_meanIntertidal_ixCrossSection = np.append(
            self.np_raster_meanIntertidal_ixCrossSection,
            point_ixCrossSection
        );
        self.np_raster_meanIntertidal_crossSection_ID = np.append(
            self.np_raster_meanIntertidal_crossSection_ID,
            self.raster_crossSection_count
        );
        
    def _append_meanHighWaterPoint_from_crossSection_to_raster(
            self, point_ixCrossSection):
        point_x = self.cc[point_ixCrossSection];
        point_y = self.rr[point_ixCrossSection];
        self.np_raster_meanHighWater_x = np.append(
            self.np_raster_meanHighWater_x,
            point_x
        );
        self.np_raster_meanHighWater_y = np.append(
            self.np_raster_meanHighWater_y,
            point_y
        );
        self.np_raster_meanHighWater_topo = np.append(
            self.np_raster_meanHighWater_topo,
            self.raster_topo_raw[point_y, point_x]
        );
        self.np_raster_meanHighWater_ixCrossSection = np.append(
            self.np_raster_meanHighWater_ixCrossSection,
            point_ixCrossSection
        );
        self.np_raster_meanHighWater_crossSection_ID = np.append(
            self.np_raster_meanHighWater_crossSection_ID,
            self.raster_crossSection_count
        );
    
    def _append_highestAstronomicalTidePoint_from_crossSection_to_raster(
            self, point_ixCrossSection):
        point_x = self.cc[point_ixCrossSection];
        point_y = self.rr[point_ixCrossSection];
        self.np_raster_highestAstronomicalTide_x = np.append(
            self.np_raster_highestAstronomicalTide_x,
            point_x
        );
        self.np_raster_highestAstronomicalTide_y = np.append(
            self.np_raster_highestAstronomicalTide_y,
            point_y
        );
        self.np_raster_highestAstronomicalTide_topo = np.append(
            self.np_raster_highestAstronomicalTide_topo,
            self.raster_topo_raw[point_y, point_x]
        );
        self.np_raster_highestAstronomicalTide_ixCrossSection = np.append(
            self.np_raster_highestAstronomicalTide_ixCrossSection,
            point_ixCrossSection
        );
        self.np_raster_highestAstronomicalTide_crossSection_ID = np.append(
            self.np_raster_highestAstronomicalTide_crossSection_ID,
            self.raster_crossSection_count
        );   

    def run_through_crossSections_over_raster(
            self,
            np_coordinatesCrossSection_refRaster_start_x,
            np_coordinatesCrossSection_refRaster_start_y,
            np_coordinatesCrossSection_refRaster_end_x,
            np_coordinatesCrossSection_refRaster_end_y):
        
        np_raster_crossSection_count = np.arange(
            0, np_coordinatesCrossSection_refRaster_start_x.size, 1,
            dtype=np.int
        );
        self._init_markerPoints_on_crossSectionProfile();
        for self.raster_crossSection_count in np_raster_crossSection_count:
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
            self.find_water_line_markers_on_crossSection(
                self.coordinatesCrossSection_refRaster_start_x,
                self.coordinatesCrossSection_refRaster_start_y,
                self.coordinatesCrossSection_refRaster_end_x,
                self.coordinatesCrossSection_refRaster_end_y
            );

    def find_water_line_markers_on_crossSection(
            self,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y):
        
        self.crossSectionProfile_topo,\
        self.rr,\
        self.cc = cs_utils.get_crossSectionProfile_from_raster(
            self.raster_topo_raw,
            coordinatesCrossSection_refRaster_start_x,
            coordinatesCrossSection_refRaster_start_y,
            coordinatesCrossSection_refRaster_end_x,
            coordinatesCrossSection_refRaster_end_y
        );
        # Low Reference Line
        lowReferenceLine_ixCrossSection = self.get_coordinatesThresholdCross(
            self.searchInterval_lowReferenceLine
        );
        if lowReferenceLine_ixCrossSection != -1:
            self._append_lowReferencePoint_from_crossSection_to_raster(
                lowReferenceLine_ixCrossSection
            );
        # Mean Low Water
        meanLowWater_ixCrossSection = self.get_coordinatesThresholdCross(
            self.searchInterval_meanLowWater
        );
        if meanLowWater_ixCrossSection != -1:
            self._append_meanLowWaterPoint_from_crossSection_to_raster(
                meanLowWater_ixCrossSection
            );
        # Mean Intertidal
        meanIntertidal_ixCrossSection = self.get_coordinatesThresholdCross(
            self.searchInterval_meanIntertidal
        );
        if meanIntertidal_ixCrossSection != -1:
            self._append_meanIntertidalPoint_from_crossSection_to_raster(
                meanIntertidal_ixCrossSection
            );
        # Mean High Water
        meanHighWater_ixCrossSection = self.get_coordinatesThresholdCross(
            self.searchInterval_meanHighWater
        );
        if meanHighWater_ixCrossSection != -1:
            self._append_meanHighWaterPoint_from_crossSection_to_raster(
                meanHighWater_ixCrossSection
            );
        # Highest Astronical Tide
        highestAstronomicalTide_ixCrossSection = self.get_coordinatesThresholdCross(
            self.searchInterval_highestAstronomicalTide
        );
        if highestAstronomicalTide_ixCrossSection != -1:
            self._append_highestAstronomicalTidePoint_from_crossSection_to_raster(
                highestAstronomicalTide_ixCrossSection
            );

    def get_coordinatesThresholdCross(self, searchInterval):
        """
        Get the coordinates of the position of the cross-section where the
         value for the topography raster is in the user-specified interval.
        
        Parameters
        ----------
        searchInterval (tuple of float)
            Interval of topography where to search for a point
            
        Returns
        -------
        numpy.int
            Index of the first detected point within the interval
            If no valid index is found, returns -1
        """
        if searchInterval is None:
            return -1;
        elif len(searchInterval) != 2:
            return -1;
        else:
            thresholdTopo_min, thresholdTopo_max = searchInterval;
            for (ixThreshold,), topo in np.ndenumerate(self.crossSectionProfile_topo):
                if topo >= thresholdTopo_min and topo <= thresholdTopo_max:
                    break;
            if ixThreshold >= (self.crossSectionProfile_topo.size-1):
                return -1;
            else:
                return ixThreshold;    
    
    def get_lowReferencePoints_over_complete_raster(self):
        np_raster_lowReferenceLine_easting = self.rEasting_origin\
            + self.np_raster_lowReferenceLine_x*self.rPixelResolution_x;
        np_raster_lowReferenceLine_northing = self.rNorthing_origin\
            + self.np_raster_lowReferenceLine_y*self.rPixelResolution_y;
        output_lowReferencePoints = np.column_stack((
            self.np_raster_lowReferenceLine_crossSection_ID,
            self.np_raster_lowReferenceLine_ixCrossSection,
            self.np_raster_lowReferenceLine_x,
            self.np_raster_lowReferenceLine_y,
            np_raster_lowReferenceLine_easting,
            np_raster_lowReferenceLine_northing,
            self.np_raster_lowReferenceLine_topo
        ));
        return pd.DataFrame(
            data=output_lowReferencePoints,
            columns=self.marker_colsHeader
        );
    
    def get_meanLowWaterPoints_over_complete_raster(self):
        np_raster_meanLowWater_easting = self.rEasting_origin\
            + self.np_raster_meanLowWater_x*self.rPixelResolution_x;
        np_raster_meanLowWater_northing = self.rNorthing_origin\
            + self.np_raster_meanLowWater_y*self.rPixelResolution_y;
        output_meanLowWaterPoints = np.column_stack((
            self.np_raster_meanLowWater_crossSection_ID,
            self.np_raster_meanLowWater_ixCrossSection,
            self.np_raster_meanLowWater_x,
            self.np_raster_meanLowWater_y,
            np_raster_meanLowWater_easting,
            np_raster_meanLowWater_northing,
            self.np_raster_meanLowWater_topo
        ));     
        return pd.DataFrame(
            data=output_meanLowWaterPoints,
            columns=self.marker_colsHeader
        );
    
    def get_meanIntertidalPoints_over_complete_raster(self):
        np_raster_meanIntertidal_easting = self.rEasting_origin\
            + self.np_raster_meanIntertidal_x*self.rPixelResolution_x;
        np_raster_meanIntertidal_northing = self.rNorthing_origin\
            + self.np_raster_meanIntertidal_y*self.rPixelResolution_y;
        output_meanIntertidalPoints = np.column_stack((
            self.np_raster_meanIntertidal_crossSection_ID,
            self.np_raster_meanIntertidal_ixCrossSection,
            self.np_raster_meanIntertidal_x,
            self.np_raster_meanIntertidal_y,
            np_raster_meanIntertidal_easting,
            np_raster_meanIntertidal_northing,
            self.np_raster_meanIntertidal_topo
        ));     
        return pd.DataFrame(
            data=output_meanIntertidalPoints,
            columns=self.marker_colsHeader
        );   
    
    def get_meanHighWaterPoints_over_complete_raster(self):
        np_raster_meanHighWater_easting = self.rEasting_origin\
            + self.np_raster_meanHighWater_x*self.rPixelResolution_x;
        np_raster_meanHighWater_northing = self.rNorthing_origin\
            + self.np_raster_meanHighWater_y*self.rPixelResolution_y;
        output_meanHighWaterPoints = np.column_stack((
            self.np_raster_meanHighWater_crossSection_ID,
            self.np_raster_meanHighWater_ixCrossSection,
            self.np_raster_meanHighWater_x,
            self.np_raster_meanHighWater_y,
            np_raster_meanHighWater_easting,
            np_raster_meanHighWater_northing,
            self.np_raster_meanHighWater_topo
        ));     
        return pd.DataFrame(
            data=output_meanHighWaterPoints,
            columns=self.marker_colsHeader
        );
    
    def get_highestAstronomicalTidePoints_over_complete_raster(self):
        np_raster_highestAstronomicalTide_easting = self.rEasting_origin\
            + self.np_raster_highestAstronomicalTide_x*self.rPixelResolution_x;
        np_raster_highestAstronomicalTide_northing = self.rNorthing_origin\
            + self.np_raster_highestAstronomicalTide_y*self.rPixelResolution_y;
        output_highestAstronomicalTidePoints = np.column_stack((
            self.np_raster_highestAstronomicalTide_crossSection_ID,
            self.np_raster_highestAstronomicalTide_ixCrossSection,
            self.np_raster_highestAstronomicalTide_x,
            self.np_raster_highestAstronomicalTide_y,
            np_raster_highestAstronomicalTide_easting,
            np_raster_highestAstronomicalTide_northing,
            self.np_raster_highestAstronomicalTide_topo
        ));     
        return pd.DataFrame(
            data=output_highestAstronomicalTidePoints,
            columns=self.marker_colsHeader
        );
    