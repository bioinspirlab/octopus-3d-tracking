%% Information
% This script is intended as an example to get you started with the
% analysis code, and reproduce some of our figures.
% See README.md for instructions, and data download links.

%% Input definitions
output_data_path = "/Users/joost/Downloads/octopus_data_output"; % this should already include the "tracking_data" subfolder
input_data_path = "/Users/joost/Downloads/octopus_data"; % Ex: /Users/joost/Downloads/octopus_data, with subfolder O15
clip = "O15_1611_15512_L1";

%% Code execution

% Get specific definitions for this clip
matpath = octo_InitializeData(clip, output_data_path, input_data_path);

% Preprocess the data - this couples the Z (depth) coordinates to the
% tracked points
preprocessdata = octo_PreProcess(clip, output_data_path, input_data_path);

octo_VisualizeCurvatureAndStrainIndividual