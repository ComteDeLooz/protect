# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt

def plot_fromPointToGroupLabel(
        raster_data, raster_shape, raster_pixel_resolution, raster_origin_refEN,
        filePath_pdPointCloud, topo_threshold_low=None, topo_threshold_high=None):
    folderPath, fileName_pdPointCloud = os.path.split(filePath_pdPointCloud);
    pdPointCloud = pd.read_csv(filePath_pdPointCloud, sep=';');
    
    easting, northing = raster_origin_refEN;
    pixel_resolution_x, pixel_resolution_y = raster_pixel_resolution;
    nRows, nColumns = raster_shape;
    extent_left = easting;
    extent_right = easting + pixel_resolution_x*nColumns;
    extent_bottom = northing + pixel_resolution_y*nRows;
    extent_top = northing;
    
    fig = plt.figure(figsize=(24,12));
    plt.imshow(
        X=raster_data,
        cmap='Greys',
        vmin=topo_threshold_low,
        vmax=topo_threshold_high,
        origin='upper',
        extent=(extent_left, extent_right, extent_bottom, extent_top)
    );
    labels_all = pdPointCloud.loc[:, 'Group_Label'].to_numpy();
    labels_unique = np.unique(labels_all);
    cmap_points = plt.cm.rainbow;
    norm = matplotlib.colors.Normalize();
    npColors = cmap_points(norm(labels_unique))
    for npIx, label in np.ndenumerate(labels_unique):
        c = npColors[npIx];
        x = pdPointCloud.loc[
            pdPointCloud['Group_Label'] == label, 'Point_Easting'
        ].to_numpy();
        y = pdPointCloud.loc[
            pdPointCloud['Group_Label'] == label, 'Point_Northing'
        ].to_numpy();
        plt.plot(x, y, marker='o', markersize=3, linestyle='', mec=c, color=c);
    plt.xlabel('Easting [m]', fontsize=12);
    plt.ylabel('Northing [m]', fontsize=12);
    mng = plt.get_current_fig_manager();
    mng.window.showMaximized();
    plt.show();
    fileName_fig = '%s.jpg' % fileName_pdPointCloud[:-4];
    fig.savefig(os.path.join(folderPath, fileName_fig), papertype='a3');
    plt.close();