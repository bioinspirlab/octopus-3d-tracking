function matpath = octo_InitializeData(clip)
% Get all the specifics that are different for clips. This program outputs
% a MAT file for easy import at the variable level in other functions

%% Set general paths
if ismac
    dropboxpath = '/Users/joost/Dropbox/MBARI/Projects/EyeRIS-octopus/';
    imgpath = '/Users/joost/Desktop/octo/exported_clips/';
    imgpath = '/Users/joost/Desktop/octo/exported_clips/ReExport_Jan25/';
%     imgpath = '/Volumes/BioInspirLab/Projects/EyeRIS/2022_EyeRIS_octopus/raytrix_exports/ReExport_19083-Jan-25/';
else
    dropboxpath = 'C:\Users\joost\Dropbox\MBARI\Projects\EyeRIS-octopus\';
    imgpath = 'C:\Users\joost\Desktop\octo\';   % Local
    imgpath = '\\atlas\BioInspirLab\Projects\EyeRIS\Octo_tracking\exported_clips\';  % Server
    imgpath = 'D:\Octo_ReExport\';  % Sabrent
    imgpath = '\\atlas\BioInspirLab\Projects\EyeRIS\Octo_tracking\ReExport_19083-Jan-25\D1457_08Octopus_20220826_161130008\'; % Server
    imgpath = 'E:\ReExport_19083-Jan-25\';
end


%% Clip-specific
guidepts = [];  % By default, all tracked points are considered accurately tracked. Guide points are those used only to constrain the spline and geodesics.
erasethese = [];
startfr = 1;
endfr = 1e4;
curvframestart = 1;
curvframeend = 1e9;
armtrimloc = 1e5;
% armtouchpoint_mm = 0;

% A note on additional parameters:
% ARMLIFT: the frame number where the arm is lifted. It is possible to
% specify a hypothetical frame number beyond the clip bounds, as long as it
% respects the frame rate. (e.g. liftoff happens 4 seconds after clip ends,
% add 4*60 = 240 to the last frame number)
% ARMTRIMLOC: The cutoff beyond which curvature data is not considered to
% determine the location of maximum curvature. Determine this by running 
% octo_VisualizeCurvatureAndStrainIndividual with mode set to normalized
% (and plotpeak and showpoints set to true), then determine in the bend
% radius plot at which y value you want the cutoff, e.g. 0.75 for 3/4 of
% the arm. This can also be defined as an absolute value of the array index
% (not recommended). Note you'll need to delete the
% segment_analysis_data_[xxx] file after changing this value to regenerate
% it.
% ARMTOUCHPOINT_MM: The location in mm from the first arm point where the
% arm last releases the substrate. Best estimated by looking in the video
% with point overlay between which tracked points this happens, and putting
% a representative value in mm using the bend radius plot from 
% octo_VisualizeCurvatureAndStrainIndividual in absolute mode.
switch clip
    case "O10_6946_L1"
        if ismac
            imgpath = '/Users/joost/Desktop/octo/exported_clips/O10_20240216/';
        else
            imgpath = 'E:\O10_20240216\';   % USB drive
%             imgpath = 'Z:\Projects\EyeRIS\2022_EyeRIS_octopus\raytrix_exports\O14_20240124\'; % atlas
        end
        imgpath = [imgpath 'D1457_03Octopus_O10_20220826_121931534' filesep];
        orthopath = [imgpath 'D1457_20220826T122519Z_Octo_ThreeDImages_6946',...
            filesep,'D1457_03Octopus_O10_20220826_121931534',...
            filesep,'D1457_20220826T121931Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T122519Z_Octo_TotalFocusOrtho_6946.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T191911Z_quad.MOV';
        quadstartframe = round(29.97*(6*60+36)+24-225/2);
        focusdepthpath = [imgpath 'D1457_20220826T122519Z_Octo_FocusedDepthMap_6946.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_O10_6946xypts.csv'];
        refpts = [7,8,9,10,11,12,13];
        ptorder = [1,2,3,4,5,6];
        startnum = 6946;
        opt.fillmode = "fill";
%         opt.fillmode = "fillmore";
        framerate = 60;
        armlift = 177;
        armtrimloc = 0.99;   % At final point
        armtouchpoint_mm = 68; % Around point 4
        gaitduration = median([41,51,69,57]);
    case "O10_7875_R1"
        if ismac
            imgpath = '/Users/joost/Desktop/octo/exported_clips/O10_20240216/';
        else
            imgpath = 'E:\O10_20240216\';   % USB drive
%             imgpath = 'Z:\Projects\EyeRIS\2022_EyeRIS_octopus\raytrix_exports\O14_20240124\'; % atlas
        end
        imgpath = [imgpath 'D1457_03Octopus_O10_20220826_121931534' filesep];
        orthopath = [imgpath 'D1457_20220826T122614Z_Octo_ThreeDImages_7875',...
            filesep,'D1457_03Octopus_O10_20220826_121931534',...
            filesep,'D1457_20220826T121931Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T122614Z_Octo_TotalFocusOrtho_7875.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T122614Z_Octo_FocusedDepthMap_7875.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T191911Z_quad.MOV';
        quadstartframe = round(29.97*(7*60+29)+15-17);
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_O10_7875xypts.csv'];
        refpts = [12,13,14];
        ptorder = [1,2,3,4,5,6,7,8,9,10,11];
        startnum = 7875;
%         numframes = 360;       % Set if you want to truncate the data
        opt.fillmode = "fill";
%         opt.fillmode = "fillmore";
        framerate = 60;
        armlift = 249;
        armtrimloc = 0.48;   % Around point 7
        armtouchpoint_mm = 82; % 2/3 of the way to point 6 from point 5
        gaitduration = median([41,51,69,57]);
    case "O14_24216_L1"
        if ismac
            imgpath = '/Users/joost/Desktop/octo/exported_clips/O14_20240124/';
        else
            imgpath = 'E:\O14_20240124\';   % USB drive
            %             imgpath = 'Z:\Projects\EyeRIS\2022_EyeRIS_octopus\raytrix_exports\O14_20240124\'; % atlas
        end
        imgpath = [imgpath 'D1457_07Octopus_O14_20220826_153322195' filesep];
        orthopath = [imgpath 'D1457_20220826T154835Z_Octo_ThreeDImages_24216',...
            filesep,'D1457_07Octopus_O14_20220826_153322195',...
            filesep,'D1457_20220826T153322Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T154835Z_Octo_TotalFocusOrtho_24216.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T223349Z_quad.MOV';
        quadstartframe = round(29.97*(15*60+3)+17);
        focusdepthpath = [imgpath 'D1457_20220826T154835Z_Octo_FocusedDepthMap_24216.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_O14_24216xypts.csv'];
        refpts = [1,2,3,4,5,6];
        ptorder = [12,35,34,11,30,10,31,9,32,8,33,7];
        guidepts = [30,31,32,33,35];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 24216;
        armlift = 2221;
        armtrimloc = 0.705;   % Halfway between point 8 and 9
        armtouchpoint_mm = 110; % Between point 9 and 10, about 1/4 of
%         the way from 9
        gaitduration = 51;  % Median arm cycle duration in seconds
    case "O14_24216_L2"
        if ismac
            imgpath = '/Users/joost/Desktop/octo/exported_clips/O14_20240124/';
        else
            imgpath = 'E:\O14_20240124\';   % USB drive
%             imgpath = 'Z:\Projects\EyeRIS\2022_EyeRIS_octopus\raytrix_exports\O14_20240124\'; % atlas
        end
        imgpath = [imgpath 'D1457_07Octopus_O14_20220826_153322195' filesep];
        orthopath = [imgpath 'D1457_20220826T154835Z_Octo_ThreeDImages_24216',...
            filesep,'D1457_07Octopus_O14_20220826_153322195',...
            filesep,'D1457_20220826T153322Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T154835Z_Octo_TotalFocusOrtho_24216.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T223349Z_quad.MOV';
        quadstartframe = round(29.97*(15*60+3)+17);
        focusdepthpath = [imgpath 'D1457_20220826T154835Z_Octo_FocusedDepthMap_24216.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_O14_24216xypts.csv'];
        refpts = [1,2,3,4,5,6];
        ptorder = [13,14,15,16,17,18,19,20,21];
        guidepts = [];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 24216;
        armlift = 6491;
        armtrimloc = 0.99;   % Location to trim output of curvature for maxima
        armtouchpoint_mm = 172;
        gaitduration = 52;
    case "O14_24216_L3"
        if ismac
            imgpath = '/Users/joost/Desktop/octo/exported_clips/O14_20240124/';
        else
            imgpath = 'E:\O14_20240124\';   % USB drive
%             imgpath = 'Z:\Projects\EyeRIS\2022_EyeRIS_octopus\raytrix_exports\O14_20240124\'; % atlas
        end
        imgpath = [imgpath 'D1457_07Octopus_O14_20220826_153322195' filesep];
        orthopath = [imgpath 'D1457_20220826T154835Z_Octo_ThreeDImages_24216',...
            filesep,'D1457_07Octopus_O14_20220826_153322195',...
            filesep,'D1457_20220826T153322Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T154835Z_Octo_TotalFocusOrtho_24216.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T223349Z_quad.MOV';
        quadstartframe = round(29.97*(15*60+3)+17);
        focusdepthpath = [imgpath 'D1457_20220826T154835Z_Octo_FocusedDepthMap_24216.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_O14_24216xypts.csv'];
        refpts = [1,2,3,4,5,6];
        ptorder = [22,27,23,28,24,25,29,26];
        guidepts = 28;
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 24216;
        armlift = 7997;
        armtrimloc = 0.786;   % Trim at point 29
        armtouchpoint_mm = 58;  % 2/3 of the way from point 24 to 25
        touchduration = (armlift-6148)/60;  % Put down approx frame 6869 from quad
        gaitduration = 44;
    case "O15_1611_5596"
        imgpath = [imgpath char(clip) filesep];
        orthopath = [imgpath 'D0_08OctopusT202208Z_EyeRIS_R26_ThreeDImages_5596',...
            filesep,'D0_08OctopusT202208Z_EyeRIS_R26_ThreeDImage_'];
        vidpath = [imgpath 'D0_08OctopusT202208Z_EyeRIS_R26_TotalFocusOrtho_5596.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_5596_armxypts.csv'];
        refpts = [11,12,13,14];
        ptorder = [1,2,3,4,5,6,7,8,9,10];
        numframes = 360;       % Set if you want to truncate the data
        opt.fillmode = "fill";
        opt.fillmode = "fillmore";
        framerate = 60;
        startfr = 1;
        endfr = 360;
    case "O15_1611_7318_L1"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T161842Z_Octo_ThreeDImages_7318',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T161842Z_Octo_TotalFocusOrtho_7318.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 9769;
        focusdepthpath = [imgpath 'D1457_20220826T161842Z_Octo_FocusedDepthMap_7318.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_7318xypts.csv'];
        refpts = [7,8,9,10,11,12,13,14,15];
        ptorder = [16,17,18,19,20,21];
        guidepts = [];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 7318;
        thresh = 3;
    case "O15_1611_7318_L2"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T161842Z_Octo_ThreeDImages_7318',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T161842Z_Octo_TotalFocusOrtho_7318.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 9769;
        focusdepthpath = [imgpath 'D1457_20220826T161842Z_Octo_FocusedDepthMap_7318.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_7318xypts.csv'];
        refpts = [7,8,9,10,11,12,13,14,15];
        ptorder = [1,2,3,4,5,6];
        guidepts = [];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 7318;

    case "O15_1611_7882_R1"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T161940Z_Octo_ThreeDImages_7882',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T161940Z_Octo_TotalFocusOrtho_7882.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 11514;
        focusdepthpath = [imgpath 'D1457_20220826T161940Z_Octo_FocusedDepthMap_7882.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_7882xypts.csv'];
        refpts = [1,2,3];
        ptorder = [16,17,27,18,19,20,21,22,23,24,25,26];
        guidepts = [];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 7882;
    case "O15_1611_7882_R3"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T161940Z_Octo_ThreeDImages_7882',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T161940Z_Octo_TotalFocusOrtho_7882.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 11514;
        focusdepthpath = [imgpath 'D1457_20220826T161940Z_Octo_FocusedDepthMap_7882.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_7882xypts.csv'];
        refpts = [1,2,3];
        ptorder = [37,4,5,6,7,8,9,14,10,15,11,13,12];
        guidepts = [13,14,15];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 7882;
        armlift = 1564;
        armtouchpoint_mm = 98;  % Near point 7
        armtrimloc = 0.6;   % Approximately at point 8
        touchduration = 13+17+(10+9)/30;  % Stationary prior to this, following step: put down ~23:20:46;20, lift off ~23:21:17;09. Note: stepped over large rock.
        gaitduration = 36.4;
    case "O15_1611_7882_R4"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T161940Z_Octo_ThreeDImages_7882',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T161940Z_Octo_TotalFocusOrtho_7882.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 11514;
        focusdepthpath = [imgpath 'D1457_20220826T161940Z_Octo_FocusedDepthMap_7882.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_7882xypts.csv'];
        refpts = [1,2,3];
        ptorder = [28,30,29,31,32,33,34,35,36];
        guidepts = [30];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 7882;
        gaitduration = 39.2;
    case "O15_1611_11589_R3"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162129Z_Octo_ThreeDImages_11589',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162129Z_Octo_TotalFocusOrtho_11589.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162129Z_Octo_FocusedDepthMap_11589.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 14805;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_11589xypts.csv'];
        refpts = [8,9,10,11,12,13];
        ptorder = [21,1,2,3,4,5,6,7];
        guidepts = [];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 11589;
        armlift = 620;
        armtouchpoint_mm = 143;
        armtrimloc = 360;   % Location to trim output of curvature for maxima
        touchduration = 27+(11)/30;  % Put down ~23:21:28;24, lift off ~ ~23:21:56;05
        gaitduration = 36.4;
    case "O15_1611_11589_R4"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162129Z_Octo_ThreeDImages_11589',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162129Z_Octo_TotalFocusOrtho_11589.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162129Z_Octo_FocusedDepthMap_11589.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 14805;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_11589xypts.csv'];
        refpts = [8,9,10,11,12,13];
        ptorder = [20,22,15,16,17,18,19];
        guidepts = [22];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 11589;
        gaitduration = 39.2;
    case "O15_1611_13209_L2"   % Also L3 and L4
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162225Z_Octo_ThreeDImages_13209',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162225Z_Octo_TotalFocusOrtho_13209.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162225Z_Octo_FocusedDepthMap_13209.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 16486;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_13209xypts.csv'];
        refpts = [5,6,7,8,9,10,11,12];
        ptorder = [1,2,3,4];
        guidepts = [];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 13209;
        gaitduration = 40.1;
    case "O15_1611_13209_L3_1"   % Stride 1
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162225Z_Octo_ThreeDImages_13209',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162225Z_Octo_TotalFocusOrtho_13209.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162225Z_Octo_FocusedDepthMap_13209.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 16486;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_13209xypts.csv'];
        refpts = [5,6,7,8,9,10,11,12];
        ptorder = [13,14,23,22,15,16,17,18,19,21,20];
        guidepts = [21,22,23];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 13209;
        armlift = 383;
        armtrimloc = 0.53;   % Halfway between 16 and 17
        armtouchpoint_mm = 66; % 1/3 of the distance from point 15 towards 16
        touchduration = 23+(26)/30;  % Put down ~23:22:24;18 (very hard to see, took middle of uncertainty window), lift off ~ ~23:22:48;14
        gaitduration = 36.4;
    case "O15_1611_13209_L3_2"   % Stride 2
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162225Z_Octo_ThreeDImages_13209',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162225Z_Octo_TotalFocusOrtho_13209.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162225Z_Octo_FocusedDepthMap_13209.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 16486;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_13209xypts.csv'];
        refpts = [5,6,7,8,9,10,11,12];
        ptorder = [13,14,23,22,15,16,17,18,19,21,20];
        guidepts = [21,22,23];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 13209;
        armlift = 2248+6*60+20*2; % Picked up 6 seconds and 20 frames
%         after frame 2248, judged from quad view
        armtrimloc = 0.76;   % Point 19
        armtouchpoint_mm = 97; % Halfway between point 17 and 18
        touchduration = (armlift-829)/60;  % Put down ~ frame 829
        gaitduration = 36.4;
    case "O15_1611_13209_L4"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162225Z_Octo_ThreeDImages_13209',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162225Z_Octo_TotalFocusOrtho_13209.mp4'];
        focusdepthpath = [imgpath 'D1457_20220826T162225Z_Octo_FocusedDepthMap_13209.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 16486;
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_13209xypts.csv'];
        refpts = [5,6,7,8,9,10,11,12];
        ptorder = [24,28,25,26,27,30,29,32,31];
        guidepts = [28,30,32];
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 13209;
        armlift = 997; 
        armtrimloc = 0.9;   % Just past point 29
        armtouchpoint_mm = 129; % 60% from 27 to 29
        touchduration = 21+(15)/30;  % Put down ~23:22:37;14 (large uncertainty), lift off ~ ~23:22:58;29
        gaitduration = 39.2;
    case "O15_1611_15512_L1"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162430Z_Octo_ThreeDImages_15512',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162430Z_Octo_TotalFocusOrtho_15512.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = 20242;
        focusdepthpath = [imgpath 'D1457_20220826T162430Z_Octo_FocusedDepthMap_15512.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_15512xypts.csv'];
        refpts = [1,9,13,14];
        ptorder = [2,3,4,5,6,7,8,10,19,17,11,12,15,18,16];
        guidepts = [8,17,19];
        opt.fillmode = "nofill";
        framerate = 60;
        startnum = 15512;
        armlift = 680;
        armtrimloc = 0.6;   % Point 10
        armtouchpoint_mm = 152; % Halfway between point 8 and point 10?
        gaitduration = 63.3;
    case "O15_1611_19083_L1"   % Used to be arm 2, then falsely labeled L2; there's a scar between 2 and 3
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162558Z_Octo_ThreeDImages_19083',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162558Z_Octo_TotalFocusOrtho_19083.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(12*60+43)+27);
        focusdepthpath = [imgpath 'D1457_20220826T162558Z_Octo_FocusedDepthMap_19083.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_19083xypts.csv'];
        refpts = [9,10,11,12,13,38];
        ptorder = [14,15,16,18,17,19];   % Arm 2 (L4)
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 19083;
        gaitduration = 63.3;
    case "O15_1611_19083_L2"   % Used to be arm1, then L3
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162558Z_Octo_ThreeDImages_19083',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162558Z_Octo_TotalFocusOrtho_19083.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(12*60+43)+27);
        focusdepthpath = [imgpath 'D1457_20220826T162558Z_Octo_FocusedDepthMap_19083.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_19083xypts.csv'];
        refpts = [9,10,11,12,13,38];
        ptorder = [1,27,2,28,43,3,29,8,39,30,4,40,34,5,31,33,6,32,7];  % Arm 1 (L3)
        guidepts = [27,28,29,30,31,32,33,34,39,40,43];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 19083;
        armlift = 950;
        armtrimloc = 225;   % Location to trim output of curvature for maxima
        armtouchpoint_mm = 58; % Approximate touch point in mm from base point
        gaitduration = 40.1;
    case "O15_1611_19083_L3"   % USed to be arm3, then L4
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162558Z_Octo_ThreeDImages_19083',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162558Z_Octo_TotalFocusOrtho_19083.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(12*60+43)+27);
        focusdepthpath = [imgpath 'D1457_20220826T162558Z_Octo_FocusedDepthMap_19083.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_19083xypts.csv'];
        refpts = [9,10,11,12,13,38];
        ptorder = [20,41,42,21,22,36,23,35,24,25,37,26];   % Arm 3 (R4)
        guidepts = [35,41,42];
        ignoreforsegment = [];  % Points to ignore for segment analysis
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 19083;
%         curvframestart = 1550;
%         curvframeend = 3200;
        armlift = 2935;
        armtrimloc = 400;   % Location to trim output of curvature for maxima
        armtouchpoint_mm = 92; % Approximate touch point in mm from base point
        touchduration = (armlift-925)/60;  % Put down at frame 925 (well visible)
        gaitduration = 36.4;
    case "O15_1611_22417_L1" % used to be called arm1, then mislabeled L2
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162728Z_Octo_ThreeDImages_22417',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162728Z_Octo_TotalFocusOrtho_22417.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(14*60+14)+18);
        focusdepthpath = [imgpath 'D1457_20220826T162728Z_Octo_FocusedDepthMap_22417.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_22417xypts.csv'];
        refpts = [8,9,10];
        ptorder = [2,1,3,4,5,6,35,7,32,33,11,34];   % Arm 1 (L2)
        guidepts = [32,33,34,35];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 22417;
    case "O15_1611_22417_L2" % L3, used to be arm2
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162728Z_Octo_ThreeDImages_22417',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162728Z_Octo_TotalFocusOrtho_22417.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(14*60+14)+18);
        focusdepthpath = [imgpath 'D1457_20220826T162728Z_Octo_FocusedDepthMap_22417.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_22417xypts.csv'];
        refpts = [8,9,10];
        ptorder = [12,13,14,15,16,17,18,19,36,23,20,37,21,22];   % Arm 2 / L3
        guidepts = [23,36,37,38];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 22417;
        armlift = 811;  
        armtrimloc = 480;   % Location to trim output of curvature for
%         maxima - around point 20
        armtouchpoint_mm = 140; % About halfway between point 18 and 19
        gaitduration = 40.1;
    case "O15_1611_22417_L3" % L4, used to be arm3
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162728Z_Octo_ThreeDImages_22417',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162728Z_Octo_TotalFocusOrtho_22417.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(14*60+14)+18);
        focusdepthpath = [imgpath 'D1457_20220826T162728Z_Octo_FocusedDepthMap_22417.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_22417xypts.csv'];
        refpts = [8,9,10];
        ptorder = [24,39,43,42,40,25,41,26,51,27,28,29,50,30,31];   % Arm 3 / L4
        guidepts = [39,41,50,51];
        ignoreforsegment = [];  % Points to ignore for segment analysis
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 22417;
        armlift = 2917;
        armtrimloc = 390;   % Location to trim output of curvature for maxima
        armtouchpoint_mm = 100;
        touchduration = (armlift-750)/60;  % Put down at frame 750
        gaitduration = 36.4;
    case "O15_1611_22417_L4"
        imgpath = [imgpath 'D1457_08Octopus_20220826_161130008' filesep];
        orthopath = [imgpath 'D1457_20220826T162728Z_Octo_ThreeDImages_22417',...
            filesep,'D1457_08Octopus_20220826_161130008',...
            filesep,'D1457_20220826T161130Z_Octo_ThreeDImage_'];
        vidpath = [imgpath 'D1457_20220826T162728Z_Octo_TotalFocusOrtho_22417.mp4'];
        quadpath = '/Users/joost/Desktop/octo/ROV/D1457_20220826T231332Z_quad.MOV';
        quadstartframe = round(29.97*(14*60+14)+18);
        focusdepthpath = [imgpath 'D1457_20220826T162728Z_Octo_FocusedDepthMap_22417.mp4'];
        datapath = [dropboxpath, 'tracking_data', filesep, 'DLTdv8_data_1611_22417xypts.csv'];
        refpts = [8,9,10];
        ptorder = [44,58,45,57,46,49,47,48,52,53,54,55,56];   % Arm 3 / L4
        guidepts = [47,49,52,54,55,57,58];
        %         numframes = 1215;       % Set if you want to truncate the data
        opt.fillmode = "fill";    % Either "nofill","fill", or "fillmore"
        framerate = 60;
        startnum = 22417;
        armlift = 1705;
        armtrimloc = 115;   % Location to trim output of curvature for maxima
        armtouchpoint_mm = 82;
        touchduration = 26+(10)/30;  % Put down ~23:27:49;02 (hard to see), lift off ~ ~23:28:15;12
        gaitduration = 39.2;
end

% Modify quadpath on MAGA
if exist("quadpath","var")
    if ispc
        quadpath = replace(quadpath,'/Users/joost/Desktop/octo/ROV/',...
            'C:\Users\joost\Desktop\octo\quad\');
        quadpath = replace(quadpath,'.MOV','.mp4');
    end
end

% A few more settings, that are identical for several arms. Here, we define
% an erasethese array. These points are included in initial depth map
% determination, but their z coordinate is later scrubbed. This is helpful
% when the depthmap is locally poorly behaved. Format: [point number, start
% frame, end frame]
switch clip
    case {"O15_1611_7882_R1","O15_1611_7882_R3","O15_1611_7882_R4"}
        erasethese = [16,500,700;...
            17,650 900;...
            27,1320,1500;...
            21,1990,2050;...
            3,280,390;...
            3,512,554;...
            24,1280,1500;...
            22,1500,1620;...
            36,1282,1430;...
            36,1823,1972;...
            ];
        thresh = 1;
    case {"O15_1611_11589_R3","O15_1611_11589_R4"}
        erasethese = [6,867,895;...
            15,1038,1220;...
            19,580,780;...
            22,500,1200;...
            ];
    case {"O15_1611_13209_L2","O15_1611_13209_L3_1","O15_1611_13209_L3_2","O15_1611_13209_L4"}
        erasethese = [1,844,900;...
            2,844,900;...
            3,844,900;...
            4,844,900;...
            15,332,396;...
            7,1,2000;...
            28,960,1300;...
            30,1200,1400;...
            26,1296,1350;...
            27,1296,1350;...
            ];
    case {"O15_1611_19083_L1","O15_1611_19083_L2","O15_1611_19083_L3"}
        erasethese = [27,2635,2720;...
            13,1854,1860;...
            13,1394,1404;...
            10,1288,1295;...
            10,1380,1385;...
            10,1430,1435;...
            9,1455,1489;...
            12,1508,1530;...
            10,1559,1565;...
            10,1737,1746;...
            10,2592,2597;...
            20,1400,1574;...
            21,1500,1727;...
            21,2400,2460;...
            26,2228,2278;...
            26,2408,2480;...
            41,1450,1640;...
            42,1450,1640;...
            4,630,660;...
            40,630,660;...
            34,630,660;...
            31,630,660;...
            6,630,660;...
            32,630,660;...
            ];
    case {"O15_1611_22417_L1","O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4"}
        erasethese = [31,550,650;...
            10,1384,1402;...
            54,500,2000;...
            55,500,2000;...
            57,500,2000;...
            58,500,2000;...
            12,1,200;...
            24,1,1200;...
            40,1,1200;...
            44,1,1200;...
            45,1,1200;...
            46,1,1200;...
            47,1,1200;...
            48,1,1200;...
            49,1,1200;...
            40,1,1200;...
            25,1,1200;...
            26,1,1200;...
            27,1,1200;...
            28,1,1200;...
            29,1,1200;...
            30,1,1200;...
            31,1,1200;...
            ];
    case {"O14_24216_L1"}
        erasethese = [12,1340,2150;...
            ];
    case {"O14_24216_L2"}
        erasethese = [13,2900,4000;...
            ];
    case {"O14_24216_L3"}
        erasethese = [18,6700,8300;...
            28,6000,8300];
end

% Add head width
clipchar = char(clip);
switch clipchar(1:3)
    case 'O10'
        headwidth = 55.6;   % mm
    case 'O14'
        headwidth = 56.6;   % mm
    case 'O15'
        headwidth = 56.9;   % mm
end

%% Save to mat file
matpath = [dropboxpath 'temp_MATLAB' filesep 'process_settings_' char(clip) '.mat'];

save(matpath);

end