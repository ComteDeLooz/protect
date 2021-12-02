# -*- coding: utf-8 -*-

import os
import numpy as np
from matplotlib import pyplot as plt

def plot_crossSectionProfile_drainageChannel(
        raster_data,
        raster_shape,
        raster_pixel_resolution,
        raster_origin_refEN,
        topo_threshold_low,
        topo_threshold_high,
        coordsCrossSection_refRaster_start,
        coordsCrossSection_refRaster_end,
        crossSection_ID,
        crossSectionProfile,
        crossSectionProfile_smooth,
        np_inflection_ixCrossSection,
        np_bottom_ixCrossSection,
        np_bottom_x,
        np_bottom_y,
        np_topLeft_ixCrossSection,
        np_topLeft_x,
        np_topLeft_y,
        np_topRight_ixCrossSection,
        np_topRight_x,
        np_topRight_y,
        folderPath_out):

    fig_map, ax_map = plt.subplots(nrows=1, ncols=1, figsize=(24,12));
    fig_cs, ax_cs = plt.subplots(nrows=1, ncols=1, figsize=(24,12));

    # Plot the raster data on the map
    easting, northing = raster_origin_refEN;
    pixel_resolution_x, pixel_resolution_y = raster_pixel_resolution;
    nRows, nColumns = raster_shape;
    extent_left = easting;
    extent_right = easting + pixel_resolution_x*nColumns;
    extent_bottom = northing + pixel_resolution_y*nRows;
    extent_top = northing;
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
    # Plotting the positions of the  points on the topography raster data
    np_bottom_easting = easting + np_bottom_x*pixel_resolution_x;
    np_bottom_northing = northing + np_bottom_y*pixel_resolution_y;
    ax_map.plot(np_bottom_easting, np_bottom_northing, 'bx', ms=3);
    np_topLeft_easting = easting + np_topLeft_x*pixel_resolution_x;
    np_topLeft_northing = northing + np_topLeft_y*pixel_resolution_y;
    ax_map.plot(np_topLeft_easting, np_topLeft_northing, 'gx', ms=3);
    np_topRight_easting = easting + np_topRight_x*pixel_resolution_x;
    np_topRight_northing = northing + np_topRight_y*pixel_resolution_y;
    ax_map.plot(np_topRight_easting, np_topRight_northing, 'gx', ms=3);
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
        ax_cs.plot(
            crossSection_x_meter,
            crossSectionProfile,
            'k',
            label='Non-smoothed topography'
        );
        ax_cs.plot(
            crossSection_x_meter,
            np.ma.masked_where(crossSectionProfile_smooth < -2, crossSectionProfile_smooth),
            'k:',
            label='Smoothed topography'
        );
    except ValueError:
        ax_cs.plot(
            crossSection_x_meter[:-1],
            crossSectionProfile,
            'k',
            label='Non-smoothed topography'
        );
        ax_cs.plot(
            crossSection_x_meter[:-1],
            np.ma.masked_where(crossSectionProfile_smooth < -2, crossSectionProfile_smooth),
            'k:',
            label='Smoothed topography'
        );
    line_style = '--';
    line_width = 1;
    for inflection_ixCrossSection in np_inflection_ixCrossSection:
        ax_cs.axvline(
            inflection_ixCrossSection * raster_pixel_resolution_mean,
            c='r', ls=':', lw=line_width,
        );
    for bottom_ixCrossSection in np_bottom_ixCrossSection:
        ax_cs.axvline(
            bottom_ixCrossSection * raster_pixel_resolution_mean,
            c='b', ls=line_style, lw=line_width,
        );
    for topLeft_ixCrossSection in np_topLeft_ixCrossSection:
        ax_cs.axvline(
            topLeft_ixCrossSection * raster_pixel_resolution_mean,
            c='g', ls=line_style, lw=line_width,
        );
    for topRight_ixCrossSection in np_topRight_ixCrossSection:
        ax_cs.axvline(
            topRight_ixCrossSection * raster_pixel_resolution_mean,
            c='g', ls=line_style, lw=line_width,
        );
    ax_cs.set_xlabel('Distance along cross-section profile [m]', fontsize=16);
    ax_cs.set_ylabel('Topography [m]', fontsize=16);
    ax_cs.set_ylim(topo_threshold_low, topo_threshold_high);
    fig_map.savefig(os.path.join(
        folderPath_out, 'CrossSection%s_map.jpg' % str(crossSection_ID)
    ));
    fig_cs.savefig(os.path.join(
        folderPath_out, 'CrossSection%s_cs.jpg' % str(crossSection_ID)
    ));
    plt.close('all');

    
    
    
    
    
    