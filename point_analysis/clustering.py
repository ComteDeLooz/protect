# -*- coding: utf-8 -*-

import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN


class FromPointToCluster(object):
    def __init__(
            self,
            pdPointCloud,
            reference_line_slope=0.4804,
            reference_line_point=(30500,204750),
            merge_clusters_interval_width_m=0.3):
        
        self.pdPointCloud = pdPointCloud;
        self.reference_line_slope = reference_line_slope;
        self.reference_line_intercept = reference_line_point;
        self.set_featureArray();

    @staticmethod
    def _distancePointToLine(
            pointEasting, pointNorthing, line_slope, line_point):
        """
        Calculate the distance from a list of points to a certain line
         (e.g. Low Reference Line or Highest Astronomical Tide)
        
        Arguments
        ---------
            pointEasting (np.array of float)
                1D list of Easting of points for which the distance to a certain line is calculated 
                No default value is set
            pointNorthing (np.array of float)
                1D list of Northing of points for which the distance to a certain line is calculated 
                No default value is set
            line_slope (float)
                Slope of the line to which the distance is calculated
                No default value is set
            line_intercept (float)
                Intercept to the northing-axis of the line to which the distance is calculated
                No default value is set
        
        Returns
        -------
            distance (np.array of float)
                1D list of distances (in meters) from the points to the line under consideration
        """
        # equation line: ax + by + c = 0
        a = line_slope;
        b = -1;
        c = -(line_slope*line_point[0] - line_point[1]);
    
        numerator = np.absolute(a*pointEasting + b*pointNorthing + c);
        denominator = np.sqrt(np.power(a,2) + np.power(b,2));
        distance = numerator / denominator;
        return distance;
    
    def set_featureArray(self):
        """Defines a numpy array with features in use to apply clustering on"""
        # Feature 1: distance to a reference line defined when instantiating class
        # Read point coordinates and transform into numpy array
        npPointCloud_easting = self.pdPointCloud.loc[:,'Point_Easting'].to_numpy();
        npPointCloud_northing = self.pdPointCloud.loc[:,'Point_Northing'].to_numpy();
        # Calculate the distance to the reference line for each point
        npDistanceToReferenceLine = self._distancePointToLine(
            npPointCloud_easting, npPointCloud_northing,
            self.reference_line_slope, self.reference_line_intercept
        );
        # Add distance to reference line information to the Pandas point cloud
        self.pdPointCloud.loc[:,'DistanceToReferenceLine'] = npDistanceToReferenceLine;
        
        # Feature 2: raw (non-smoothed) topography at the point
        npTopo = self.pdPointCloud.loc[:,'Point_Topography'].to_numpy();
        
        # Define the feature array to be used in the clustering
        self.npFeature = np.column_stack((npDistanceToReferenceLine, npTopo));

    def apply_clustering_dbscan(self, eps=0.10, min_samples=35):
        self.pdPointCloud['Cluster_Label'] = -1;
    
        scaler = StandardScaler();
        self.npFeature_scaled = scaler.fit_transform(self.npFeature);
        
        clustering = DBSCAN(eps, min_samples).fit(self.npFeature_scaled);
        self.pdPointCloud.loc[:,'Cluster_Label'] = clustering.labels_;
        self.pdPointCloud.loc[:,'Cluster_Is_Core_Sample'] = 0;
        cluster_core_samples = clustering.core_sample_indices_.astype(np.int);
        ix = self.pdPointCloud.iloc[cluster_core_samples,:].index;
        self.pdPointCloud.loc[ix, 'Cluster_Is_Core_Sample'] = 1;
    
    def merge_clusters_with_comparable_topography(self, interval_width_m=0.3):
        """
        The clustering might result in a division of points in two separate
            clusters that are better to group in one.
        This is often the case when topography is quite similar but the distance
            to the reference line is more different.
        This function groups two clusters in one if the median topography between
            both is not too different
        
        Parameter
        ---------
        interval_width_m (float), optional
            Maximum distance in median topography for clusters to be grouped in one
        """
        # Get the number of clusters in the point cloud (excluding -1)
        num_clusters = self.pdPointCloud.loc[:,'Cluster_Label'].max() + 1;
        # Create a numpy array with the unique cluster labels
        np_clusters_label = np.arange(num_clusters);
        # Initialize a numpy array storing the median topography of the clusters
        np_clusters_median_topo = np.zeros((num_clusters,), dtype=np.float32);
        # Fill in the median topography of the clusters in the numpy array
        for cluster_label in np_clusters_label:
            cluster_median_topo = self.pdPointCloud.loc[
                self.pdPointCloud['Cluster_Label'] == cluster_label,
                'Point_Topography'
            ].median();
            np_clusters_median_topo[cluster_label] = cluster_median_topo;
        # Loop through the numpy array with median topography per cluster
        # Initially, the index of the element in the array marks the original
        #  cluster label
        for (npIx,), cluster_median_topo in np.ndenumerate(np_clusters_median_topo):
            # If the last cluster label is encountered, no elements remain to
            #  compare with
            if npIx == (np_clusters_median_topo.size - 1):
                break;
            # If the index of the current element of the median topography array
            #  does not correspond with the cluster label, the cluster label is
            #  already updated in a previous run through the loop. Therefore, this
            #  element is not reprocessed
            if npIx != np_clusters_label[npIx]:
                continue;
            # Make a hard copy of the median topography array
            np_clusters_median_topo_check = np.copy(np_clusters_median_topo);
            # Calculate the absolute difference in median topography between
            #  each cluster element and the current element
            np_clusters_median_topo_diff = np.absolute(
                np_clusters_median_topo_check - cluster_median_topo
            );
            # Check for each difference if it falls within the interval specified
            #  as parameter to this method
            np_clusters_median_topo_comp = np_clusters_median_topo_diff\
                                           <= interval_width_m
            # Find the indices of the array where the median topography difference
            #  is within the interval. Only the indices after the current index
            #  element are taken into account because the previous elements were
            #  already compared in a prevous run of the loop
            npIx_merge_candidates = np.argwhere(
                np_clusters_median_topo_comp[(npIx + 1):]
            )[:,0];
            # If an index is found where the median topography is within the interval,
            #  replace the cluster element label with the current label 
            if not npIx_merge_candidates.size == 0:
                np_clusters_label[(npIx_merge_candidates + npIx + 1)] = npIx;
        # Store the information in the Pandas DataFrame in a separate column
        for (old_label,), new_label in np.ndenumerate(np_clusters_label):
            self.pdPointCloud.loc[
                self.pdPointCloud['Cluster_Label'] == old_label,
                'Cluster_Label_Merge'
            ] = new_label;
        self.pdPointCloud.loc[
            self.pdPointCloud['Cluster_Label'] == -1,
            'Cluster_Label_Merge'
            ] = -1;
    
    def reorganize_clusters_found(
            self, length_beach_area_pxl, spacing_between_profiles_pxl,
            min_spatial_extent_pct=60):
        """
        Reorganise the recorded cluster order to get the same order for all
            raster data over time
        """
        # Get the number of clusters in the point cloud (excluding -1)
        try:
            cluster_col_name = 'Cluster_Label_Merge';
            np_cluster_labels_old = self.pdPointCloud.loc[
                :, cluster_col_name
            ].unique();
        except KeyError: # If the column Cluster_Label_Merge does not exist
            cluster_col_name = 'Cluster_Label';
            np_cluster_labels_old = self.pdPointCloud.loc[
                :, cluster_col_name
            ].unique();
        # Exclude the points not allocated to a cluster (-1)
        np_cluster_labels_old = np_cluster_labels_old[
            np_cluster_labels_old != -1
        ];
        np_cluster_labels_old = np_cluster_labels_old.astype(np.uint);
        # Count the number of unique cluster labels
        num_clusters = np_cluster_labels_old.size;
        np_cluster_info = np.vstack((
            np_cluster_labels_old,
            np.zeros((num_clusters,), dtype=np.uint),
            np.zeros((num_clusters,), dtype=np.float32)
        ));
        # Loop over the different original clusters
        for cluster_info in np_cluster_info.T:
            # Filter the Pandas DataFrame to keep only the current cluster points
            pdPointCloud_perCluster = self.pdPointCloud.loc[
                self.pdPointCloud[cluster_col_name] == cluster_info[0], :
            ];
            # Count the number of points in the current cluster
            cluster_info[1] = pdPointCloud_perCluster.shape[0];
            # Calculate the median of the topography of all points within the
            #  current cluster
            cluster_info[2] = pdPointCloud_perCluster.loc[
                :, 'Point_Topography'
            ].median();
        # Sort the 2D array from low to high median topography
        np_cluster_info = np_cluster_info[
            :, np_cluster_info[2,:].argsort()
        ];
        label_reorganize = 1;
        self.pdPointCloud['Cluster_Label_Reorganize'] = -1;
        for (ixCol,), old_label in np.ndenumerate(np_cluster_info[0,:]):
            if np_cluster_info[1, ixCol]\
                    <= min_spatial_extent_pct/100 * length_beach_area_pxl/spacing_between_profiles_pxl:
                continue;
            self.pdPointCloud.loc[
                self.pdPointCloud[cluster_col_name] == old_label,
                'Cluster_Label_Reorganize'
            ] = label_reorganize * 100;
            label_reorganize += 1;
    
    def get_point_cloud_with_appended_information(self):
        return self.pdPointCloud;
    