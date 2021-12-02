# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_label_colour(label):
    if label == 100:
        colour = 'r';
    elif label == 200:
        colour = 'g';
    elif label == 300:
        colour = 'b';
    elif label == 400:
        colour = 'c';
    elif label == 500:
        colour = 'y';
    elif label == 600:
        colour = 'm';
    elif label == -1:
        colour = 'k';
    return colour;

def plot_fromPointToCluster(
        raster_data, raster_shape, raster_pixel_resolution, raster_origin_refEN,
        filePath_pdPointCloud, topo_threshold_low=None, topo_threshold_high=None,
        reference_line_slope=0, reference_line_point=(305000,204750)):
    fig_map_raster_only, ax_map_raster_only = plt.subplots(
        nrows=1, ncols=1, figsize=(24,12)
    );
    fig_map, ax_map = plt.subplots(
        nrows=1, ncols=1, figsize=(24,12)
    );
    fig_dbscan, ax_dbscan = plt.subplots(
        nrows=1, ncols=1, figsize=(24,12)
    );
    
    folderPath, fileName_pdPointCloud = os.path.split(filePath_pdPointCloud);
    pdPointCloud = pd.read_csv(filePath_pdPointCloud, sep=';');
    
    # Plot the raster data on the map
    easting, northing = raster_origin_refEN;
    pixel_resolution_x, pixel_resolution_y = raster_pixel_resolution;
    nRows, nColumns = raster_shape;
    extent_left = easting;
    extent_right = easting + pixel_resolution_x*nColumns;
    extent_bottom = northing + pixel_resolution_y*nRows;
    extent_top = northing;
    im_raster_only = ax_map_raster_only.imshow(
        X=raster_data,
        cmap='Spectral',
        vmin=topo_threshold_low,
        vmax=topo_threshold_high,
        origin='upper',
        extent=(extent_left, extent_right, extent_bottom, extent_top))
    ax_map.imshow(
        X=raster_data,
        cmap='Greys',
        vmin=topo_threshold_low,
        vmax=topo_threshold_high,
        origin='upper',
        extent=(extent_left, extent_right, extent_bottom, extent_top)
    );
    
    # Plot the reference line on the map
    a = reference_line_slope;
    b = -1;
    c = -(
        reference_line_slope*reference_line_point[0] - reference_line_point[1]
    );
    line_easting_start = extent_left;
    line_northing_start = (-a*line_easting_start - c) / b;
    line_northing_end = extent_top;
    line_easting_end = (-b*line_northing_end - c)/a;
    ax_map.plot(
        [line_easting_start, line_easting_end],
        [line_northing_start, line_northing_end],
        'r--'
    );
    
    clusters_all = pdPointCloud.loc[:,'Cluster_Label_Reorganize'].to_numpy();
    clusters_unique = np.unique(clusters_all);
    for npIx, cluster in np.ndenumerate(clusters_unique):
        c = get_label_colour(cluster)
        x = pdPointCloud.loc[
            pdPointCloud['Cluster_Label_Reorganize'] == cluster,
            'Point_Easting'
        ].to_numpy();
        y = pdPointCloud.loc[
            pdPointCloud['Cluster_Label_Reorganize'] == cluster, 'Point_Northing'
        ].to_numpy();
        # Plot the points on the map within the color of the cluster
        ax_map.plot(
            x, y, marker='o', markersize=1, linestyle='', mec=c, color=c
        );
        # Plot the points on the feature graph with the same color as on the map
        dist_cluster_only_core_samples = pdPointCloud.loc[
            (pdPointCloud['Cluster_Label_Reorganize'] == cluster)\
            & (pdPointCloud['Cluster_Is_Core_Sample'] == 1),
            'DistanceToReferenceLine'
        ].to_numpy();
        dist_cluster_all = pdPointCloud.loc[
            pdPointCloud['Cluster_Label_Reorganize'] == cluster,
            'DistanceToReferenceLine'
        ].to_numpy();
        topo_cluster_only_core_samples = pdPointCloud.loc[
            (pdPointCloud['Cluster_Label_Reorganize'] == cluster)\
            & (pdPointCloud['Cluster_Is_Core_Sample'] == 1),
            'Point_Topography'
        ].to_numpy();
        topo_cluster_all = pdPointCloud.loc[
            pdPointCloud['Cluster_Label_Reorganize'] == cluster,
            'Point_Topography'
        ].to_numpy();
        ax_dbscan.plot(
            dist_cluster_all, topo_cluster_all,
            marker='o', markersize=3, linestyle='', mec=c, color=c
        );
        ax_dbscan.plot(
            dist_cluster_only_core_samples, topo_cluster_only_core_samples,
            marker='o', markersize=3, linestyle='', mec='k', color=c
        );
        if cluster != -1:
            ax_dbscan.axvline(
                x=np.median(dist_cluster_all), c=c, ls='--'
            );
            ax_dbscan.axvline(
                x=np.median(dist_cluster_all) + np.std(dist_cluster_all),
                c=c, ls=':'
            );
            ax_dbscan.axvline(
                x=np.median(dist_cluster_all) - np.std(dist_cluster_all),
                c=c,ls=':'
            );
            ax_dbscan.axhline(
                y=np.median(topo_cluster_all), c=c, ls='--'
            );
            ax_dbscan.axhline(
                y=np.median(topo_cluster_all) + np.std(topo_cluster_all),
                c=c, ls=':'
            );
            ax_dbscan.axhline(
                y=np.median(topo_cluster_all) - np.std(topo_cluster_all),
                c=c, ls=':'
            );
    divider = make_axes_locatable(ax_map_raster_only);
    cax = divider.append_axes('right', size='5%', pad=0.05);
    cbar = plt.colorbar(im_raster_only, cax=cax);
    cbar.set_label('Topography [m]', fontsize=16);
    ax_map_raster_only.tick_params(axis='x', labelsize=16, labelrotation=90);
    ax_map_raster_only.tick_params(axis='y', labelsize=16);
    ax_map_raster_only.set_xlabel('Easting [m]', fontsize=16);
    ax_map_raster_only.set_ylabel('Northing [m]', fontsize=16);
    fileName_fig_map_raster_only = '%s_map_raster.jpg' % fileName_pdPointCloud[:-4];
    fig_map_raster_only.savefig(
        os.path.join(folderPath, fileName_fig_map_raster_only), papertype='a4'
    );
    
    ax_map.tick_params(axis='x', labelsize=16, labelrotation=90);
    ax_map.tick_params(axis='y', labelsize=16);
    ax_map.set_xlabel('Easting [m]', fontsize=16);
    ax_map.set_ylabel('Northing [m]', fontsize=16);
    fileName_fig_map = '%s_map.jpg' % fileName_pdPointCloud[:-4];
    fig_map.savefig(
        os.path.join(folderPath, fileName_fig_map), papertype='a4'
    );
    
    ax_dbscan.set_xlabel('Distance To Reference Line [m]', fontsize=12);
    ax_dbscan.set_ylabel('Point Topography [m]', fontsize=12);
    fileName_fig_dbscan = '%s_dbscan.jpg' % fileName_pdPointCloud[:-4];
    fig_dbscan.savefig(
        os.path.join(folderPath, fileName_fig_dbscan), papertype='a4'
    );
    plt.close('all');
    