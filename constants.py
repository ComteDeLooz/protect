# -*- coding: utf-8 -*-

# Filters out "side-effects" for the beach area : int
BEACH_FILTERING_BINARY_CLOSING_ITERATIONS = 100;

# Constants related to the interdal bar analysis
BAR = {
    # Length (in meter) of the median filter when searching for inflection points : int 
    'INFLECTION_POINTS_MEDIAN_FILTER_KERNEL_METER': 11
};

# Constants related to the drainage channel analysis
CHANNEL = {
    # Length (in meter) of the median filter when searching for inflection points : int
    'INFLECTION_POINTS_MEDIAN_FILTER_KERNEL_METER': 11
};

# Minimum number of feature points on the raster to continue the analysis : int
MINIMUM_NUMBER_POINTS_ON_RASTER_TO_DO_ANALYSIS = 50;
