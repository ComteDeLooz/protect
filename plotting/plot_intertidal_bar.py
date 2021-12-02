# -*- coding: utf-8 -*-

import os
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict

def plot_crossSectionProfile_intertidalBar(
        raster_data,
        raster_shape,
        raster_pixel_resolution,
        raster_origin_refEN,
        topo_interval,
        coordsCrossSection_refRaster_start,
        coordsCrossSection_refRaster_end,
        crossSection_ID,
        crossSectionProfile,
        crossSectionProfile_smooth,
        np_inflection_pointID,
        np_inflection_ixCrossSection,
        np_crestTrough_pointID,
        np_crest_ixCrossSection,
        np_crest_x,
        np_crest_y,
        np_trough_ixCrossSection,
        np_trough_x,
        np_trough_y,
        folderPath_out):
    
    # Initialize figures
    fig_map, ax_map = plt.subplots(nrows=1, ncols=1, figsize=(24,12));
    fig_cross, ax_cross = plt.subplots(nrows=1, ncols=1, figsize=(24,12));
    fig_cross.suptitle('Cross-Section Profile %s' % crossSection_ID);
    # Plotting topography raster data
    easting, northing = raster_origin_refEN;
    pixel_resolution_x, pixel_resolution_y = raster_pixel_resolution;
    nRows, nColumns = raster_shape;
    extent_left = easting;
    extent_right = easting + pixel_resolution_x*nColumns;
    extent_bottom = northing + pixel_resolution_y*nRows;
    extent_top = northing;
    topo_threshold_low, topo_threshold_high = topo_interval;
    ax_map.imshow(
        X=raster_data,
        cmap='Greys',
        vmin=topo_threshold_low,
        vmax=topo_threshold_high,
        origin='upper',
        extent=(extent_left, extent_right, extent_bottom, extent_top)
    );
    # Plotting cross-section profile line on the topography raster data
    coordsCrossSection_refRaster_start_x,\
    coordsCrossSection_refRaster_start_y = coordsCrossSection_refRaster_start;
    coordsCrossSection_refRaster_end_x,\
    coordsCrossSection_refRaster_end_y = coordsCrossSection_refRaster_end;
    coordsCrossSection_start_easting = easting\
                                       + coordsCrossSection_refRaster_start_x*pixel_resolution_x;
    coordsCrossSection_start_northing = northing\
                                        + coordsCrossSection_refRaster_start_y*pixel_resolution_y;
    coordsCrossSection_end_easting = easting\
                                     + coordsCrossSection_refRaster_end_x*pixel_resolution_x;
    coordsCrossSection_end_northing = northing\
                                      + coordsCrossSection_refRaster_end_y*pixel_resolution_y;
    ax_map.plot(
        [coordsCrossSection_start_easting, coordsCrossSection_end_easting],
        [coordsCrossSection_start_northing, coordsCrossSection_end_northing],
        'k:'
    );
    # Plotting the positions of the crest and trough points on the topography raster data
    np_crest_easting = easting + np_crest_x*pixel_resolution_x;
    np_crest_northing = northing + np_crest_y*pixel_resolution_y;
    ax_map.plot(np_crest_easting, np_crest_northing, 'gx', ms=3);
    np_trough_easting = easting + np_trough_x*pixel_resolution_x;
    np_trough_northing = northing + np_trough_y*pixel_resolution_y;
    ax_map.plot(np_trough_easting, np_trough_northing, 'bx', ms=3);
    # Naming the axes of the subplot
    ax_map.set_xlabel('Easting [m]', fontsize=12);
    ax_map.set_ylabel('Northing [m]', fontsize=12);
    
    # Transform the index on profile into a distance measure (meter)
    raster_pixel_resolution_mean = np.mean(np.absolute(
        [pixel_resolution_x, pixel_resolution_y]
    ));
    crossSection_x_meter = np.arange(
        0,
        crossSectionProfile.size * raster_pixel_resolution_mean,
        raster_pixel_resolution_mean
    );
    try:
        ax_cross.plot(
            crossSection_x_meter,
            crossSectionProfile,
            'k',
            label='Non-smoothed topography'
        );
        ax_cross.plot(
            crossSection_x_meter,
            np.ma.masked_where(crossSectionProfile_smooth < -2, crossSectionProfile_smooth),
            'k:',
            label='Smoothed topography'
        );
    except ValueError:
        ax_cross.plot(
            crossSection_x_meter[:-1],
            crossSectionProfile,
            'k',
            label='Non-smoothed topography'
        );
        ax_cross.plot(
            crossSection_x_meter[:-1],
            np.ma.masked_where(crossSectionProfile_smooth < -2, crossSectionProfile_smooth),
            'k:',
            label='Smoothed topography'
        );
    # Put vertical lines along the cross-section profile where crest and trough
    #  points are found (both strict conditions and not)
    for (point_count,), crestTrough_pointID in np.ndenumerate(np_crestTrough_pointID):
        inflection_ixCrossSection = np_inflection_ixCrossSection[
            np_inflection_pointID==crestTrough_pointID
        ];
        crest_ixCrossSection = np_crest_ixCrossSection[point_count];
        trough_ixCrossSection = np_trough_ixCrossSection[point_count];
        # If inflection point falls together with trough point, subtract 1 to
        #  make a distinction on the figure
        if inflection_ixCrossSection == trough_ixCrossSection:
            inflection_ixCrossSection -= 1;
        elif inflection_ixCrossSection == crest_ixCrossSection:
            inflection_ixCrossSection += 1;
        line_style = '--';
        line_width = 1;
        ax_cross.axvline(
            inflection_ixCrossSection * raster_pixel_resolution_mean,
            c='r', ls=':', lw=line_width,
            label='Inflection'
        );
        ax_cross.axvline(
            crest_ixCrossSection * raster_pixel_resolution_mean,
            c='g', ls=line_style, lw=line_width,
            label='Crest'
        );
        ax_cross.axvline(
            trough_ixCrossSection * raster_pixel_resolution_mean,
            c='b', ls=line_style, lw=line_width,
            label='Trough'
        );
        ax_cross.set_xlabel('Distance along cross-section profile [m]', fontsize=16);
        ax_cross.set_ylabel('Topography [m]', fontsize=16);
        ax_cross.set_ylim(topo_threshold_low, topo_threshold_high);
        handles, labels = ax_cross.get_legend_handles_labels();
        by_label = OrderedDict(zip(labels, handles));
        ax_cross.legend(by_label.values(), by_label.keys(), fontsize=16);
    fig_map.savefig(os.path.join(
        folderPath_out, 'CrossSection%s_map.jpg' % str(crossSection_ID)
    ));
    fig_cross.savefig(os.path.join(
        folderPath_out, 'CrossSection%s.jpg' % str(crossSection_ID)
    ));
    plt.close('all');