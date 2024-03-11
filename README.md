# octopus-3d-tracking
MATLAB code for analyzing 3D image data from the deep-sea plenoptic imaging instrument [EyeRIS](https://www.mbari.org/technology/eyeris/).

This code was used to analyze the data presented in a manuscript that's currently under review. Data needed to run the code can be found on Zenodo, under DOI [10.5281/zenodo.10795493](https://doi.org/10.5281/zenodo.10795493).

## Overview
These analysis programs assume you have depth maps in tiff format, as well as tracked points in CSV format, as generated by the [DLTdv8 MATLAB app](https://biomech.web.unc.edu/dltdv/) by Tyson Hedrick (DOI [10.1088/1748-3182/3/3/034001](https://doi.org/10.1088/1748-3182/3/3/034001)).
Important definitions related to the datasets (referred to as 'clips') are included in "octo_InitializeData.m". The basic process is to do the heavy processing, i.e. the reading, filtering, and applying depth maps, as well as data smoothing, prior to visualization and data exports. This all happens in "octo_PreProcess.m", which is called by most of the analysis and export programs - it will attempt to cache the preprocessed files - be sure to clear this cache folder ("temp_MATLAB" by default) if you've made changes to the point data and/or octo_InitializeData file.

## Contact

## Instructions for use

* Download the data you'd like to use.
* Run "octo_InitializeData.m"; it will create a definitions file "octo_DataLocations.m".
* Enter the appropriate paths in "octo_DataLocations.m".
* Optional: preprocess the clips of interest at a convenient time. This can take a while. "octo_batch.m" is a convenient way to batch process multiple clips.
* Run a visualization, analysis or export script of your choosing. Most are implemented as functions, and require the clip name as defined in "octo_InitializeData.m" (e.g., "O15_1611_19083_L3").

## Additional details

### Analysis programs

* octo_PreProcess.m: Used to calculate 3D point data from depth maps and 2D point tracking. Results are cached in the cache folder.
* octo_AnalyzeSegments.m: Calculates segment properties. The "curvPeaks" function inside this file optionally calculates the location of maximum curvature - this requires certaina dditional settings in "octo_InitializeData.m"
* calcGeodesic.m: An attempt at tracking the arms across the depth map surface, rather than interpolation depth between points with a spline. Not currently used.
* calcPath.m: Similar to calcGeodesic.m
* getSphereFrom3Points.m: Used to calculate bend radius.
* octo_ProcessDepthMap.m: Used to process depth maps prior to determining associated z coordinates from 2D data.
* surfaceSpline.m: Determines spline interpolation between tracked points.

### Visualizations

* octo_VisualizeCurvatureAndStrainIndividual.m: The go-to program to generate strain and curvature plots for a single clip.
* octo_VisualizeCurvatureAggregate.m: Aggregate statistics for curvature across selected clips. 
* octo_VisualizeGlobalMovement.m: Read stride observations from a Google sheet and generates a visualization.
* octo_VisualizePeakCurvatureAggregate.m: Peak curvature across selected clips. 
* octo_VisualizeStrainAggregate.m: Aggregate statistics for strain across selected clips. 
* octo_DataVisualization_Other.m: Various little visualization experiments in one program.

### Video export programs

* octo_ExportVideo_3Dviews.m: Generate a video with two different angles of the spline-interpolated arm, with total focus imagery on the background.
* octo_ExportVideo_OverlayPoints.m: Overlay points on videos.
* octo_ExportVideo_OverlayPointsAndSpline.m: Generate a video of total focus and depth map, with points and spline overlaid.
* octo_ExportVideo_Quad3D.m: Generate a multi-panel video with two ROV views of the animal, a depth map from EyeRIS data, total focus image with points and spline overlaid, and a 3D plot of the interpolated arm attitude.
* octo_ExportVideoPNG_DepthMap.m: Export depth map images as PNGs.
* octo_ExportVideoPNG_TotalFocus_OverlayPoints.m: Export total focus images with point overlaid (PNG and/or MP4). Can include multiple arms for the same sequence.

### Misc programs

* convert4DtiffsTo1D.m: Default export processes of the rays files generate 4D tiffs. However, x/Y values can be inferred from just a single reference tiff. This script copies z (depth) coordinate of the TIFFs over to a new folder for reduced size and no effective loss of information.