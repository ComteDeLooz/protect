# -*- coding: utf-8 -*-

import numpy as np
from skimage import draw

def get_coordinates_crossSections_refRaster(
        raster_pixel_resolution,
        raster_shape,
        profile_orientation,
        profile_spacing_m=2):
    """
    Calculates the start and end coordinates of the cross-sections (profiles)
    oriented at a certain angle with the x-axis and expressed in image coordinates.
    The orientation of the axes is according the image coordinates system
    (x: + in horizontal right direction; y: + in vertical down direction).
    The origin of the coordinate system is at the top left of the image.
    
    Parameters
    ----------    
    raster_pixel_resolution (tuple of float)
        (pixel_resolution_x, pixel_resolution_y)
    raster_shape (tuple of int)
        (number_of_rows, number_of_columns)
    profile_orientation (float)
        Absolute angle of the profile to be selected from the raster,
        expressed in degrees:
            = 0: corresponds to the x-axis (horizontal profile)
            = 90: corresponds to the y-axis (vertical profile)
            + orientation: clock-wise
            - orientation: counter clock-wise
    profile_spacing_m (float), optional
        Perpendicular spacing [meter] between subsequent profile lines.
        Default is 2 meter.

    Returns
    -------
    profile_start_x (numpy.array of int)
        All subsequent starting coordinates (x-axis) of the profiles
    profile_start_y (numpy.array of int)
        All subsequent starting coordinates (y-axis) of the profiles 
    profile_end_x (numpy.array of int)
        All subsequent ending coordinates (x-axis) of the profiles 
    profile_end_y (numpy.array of int)
        All subsequent ending coordinates (y-axis) of the profiles
    """
    rPixelResolution_x, rPixelResolution_y = raster_pixel_resolution;
    rPixelResolution_mean = np.mean(
        np.absolute([rPixelResolution_x, rPixelResolution_y])
    );
    profileSpacing_pxl = np.int(
        np.around(profile_spacing_m / rPixelResolution_mean)
    );
    
    raster_nRows, raster_nCols = raster_shape;
    raster_x_max = raster_nCols - 1; 
    raster_y_max = raster_nRows - 1;
    
    # Profiles correspond to rows of the raster (horizontal lines)
    if profile_orientation == 0:
        # np.array with all x-coordinates of the starting coordinate of the profiles 
        profile_start_x = np.full(
            (np.floor_divide(raster_nRows, profileSpacing_pxl),),
            0, dtype=np.int
        );
        # np.array with all x-coordinates of the ending coordinate of the profiles
        profile_end_x = np.full(
            (np.floor_divide(raster_nRows, profileSpacing_pxl),),
            (raster_nCols-1), dtype=np.int
        );
        # np.array with all y-coordinates of the starting coordinate of the profiles 
        profile_start_y = np.arange(
            0, raster_nRows, profileSpacing_pxl, dtype=np.int
        );
        # np.array with all y-coordinates of the ending coordinate of the profiles 
        profile_end_y = np.arange(
            0, raster_nRows, profileSpacing_pxl, dtype=np.int
        );
    
    # Profiles correspond to columns of the raster (vertical lines)
    elif profile_orientation == 90:
        profile_start_y = np.full(
            (np.floor_divide(raster_nCols, profileSpacing_pxl),),
            0, dtype=np.int
        );
        profile_end_y = np.full(
            (np.floor_divide(raster_nCols, profileSpacing_pxl),),
            (raster_nRows-1), dtype=np.int
        );
        profile_start_x = np.arange(
            0, raster_nCols, profileSpacing_pxl, dtype=np.int
        );
        profile_end_x = np.arange(
            0, raster_nCols, profileSpacing_pxl, dtype=np.int
        );
    
    else:
        profile_start_x = np.array([], dtype=np.int);
        profile_start_y = np.array([], dtype=np.int);
        profile_end_x = np.array([], dtype=np.int);
        profile_end_y = np.array([], dtype=np.int);
        
        # y = ax + b
        rico = np.tan(np.radians(profile_orientation));
        
        # Determine the minimum and maximum offset (y-axis) of the straight line
        #   with the defined orientation. 
        if profile_orientation > 0:
            # y = ax + b (x=0; y=max); origin upper left
            b_max = raster_y_max;
            # y = ax + b (x=max; y=0); origin upper left
            b_min = - rico*raster_x_max;
        else:
            # y = ax + b (x=max; y=max); origin upper left
            b_max = raster_y_max - rico*raster_x_max;
            # y = ax + b (x=0; y=0); origin upper left
            b_min = 0;
        # Calculating the spacing in offset between subsequent profile lines
        #   based on perpendicular spacing provided as input.
        offsetSpacing = np.int(
            profileSpacing_pxl / np.cos(np.radians(profile_orientation))
        );
        if offsetSpacing < 1:
            offsetSpacing = 1;
        # Loop through all offsets, starting with maximum
        #   A new straight line is calculated for each pixel in vertical (y) direction
        for b in np.arange(b_max, b_min, -offsetSpacing):
            if profile_orientation > 0:
                # y = ax + b (y=0); origin upper left
                start_x = -b/rico;
                # y = ax + b (x=0); origin upper left
                start_y = b;
                
                if start_x < 0:
                    start_x = 0;
                if start_y < 0:
                    start_y = 0;  
                
                # y = ax + b (y=max); origin upper left
                end_x = (raster_y_max-b) / rico;
                # y = ax + b (x=max); origin upper left
                end_y = rico*raster_x_max + b;
                
                if end_x > raster_x_max:
                    end_x = raster_x_max;
                if end_y > raster_y_max:
                    end_y = raster_y_max;
            else:
                # y = ax + b (y=max); origin upper left
                start_x = (raster_y_max - b) / rico; 
                # y = ax + b (x=0); origin upper left
                start_y = b;
                
                # The minimum possible value for the x-coordinate is zero.
                if start_x < 0:
                    start_x = 0;
                # The maximum possible value for the y-coordinate is (number of rows - 1)
                if start_y > raster_y_max:
                    start_y = raster_y_max;
                
                # y = ax + b (y=0); origin upper left
                end_x = - b/rico;
                # y = ax + b (x=max); origin upper left
                end_y = rico*raster_x_max + b;
            
                if end_x > raster_x_max:
                    end_x = raster_x_max;
                if end_y < 0:
                    end_y = 0;
            
            profile_start_x = np.append(profile_start_x, np.int(start_x));
            profile_start_y = np.append(profile_start_y, np.int(start_y));
            profile_end_x = np.append(profile_end_x, np.int(end_x));
            profile_end_y = np.append(profile_end_y, np.int(end_y));
                
    return profile_start_x, profile_start_y, profile_end_x, profile_end_y;

def get_crossSectionProfile_from_raster(
        raster,
        coordinatesCrossSection_refRaster_start_x,
        coordinatesCrossSection_refRaster_start_y,
        coordinatesCrossSection_refRaster_end_x,
        coordinatesCrossSection_refRaster_end_y):
    """
    Calculates the coordinates of the pixels (row, column) that belong to
    the line corresponding with the provided cross-section. In addition,
    the topography values over the cross-section are also extracted.
    """
    # Indices of pixels that belong to the line.
    rr, cc = draw.line(
        coordinatesCrossSection_refRaster_start_y,
        coordinatesCrossSection_refRaster_start_x,
        coordinatesCrossSection_refRaster_end_y,
        coordinatesCrossSection_refRaster_end_x
    );
    crossRaster = raster[rr, cc];
    return crossRaster, rr, cc;
