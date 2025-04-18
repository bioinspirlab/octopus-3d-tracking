% Preprocess depth map
function [depthmap,opt] = octo_ProcessDepthMap(i,orthopath,startnum,opt)

if nargin<4
    opt.fillmode="nofill";
end

if ~isfield(opt,"med3")
    opt.med3 =  1; % Number of frames to use (median filter) Default 5
end

% Definitions
med3 = opt.med3;  
med2 = 11;  % 2D median filter window size

% Determine images to load
imnum = i-(med3-1)/2 : i+(med3-1)/2;
imnum(imnum<1) = [];

% Load from cache if possible
if isfield(opt,"cache")
    [a,b] = ismember(imnum,opt.cachenum);
    sz = size(opt.cache);
    depthmap = nan(sz(1),sz(2),numel(imnum));
    depthmap(:,:,a) = opt.cache(:,:,b(b>0));
    todo = find(~a);
else
    todo = 1:numel(imnum);
end


for k = todo
    try
        imname_depth = [orthopath,num2str(startnum-1+imnum(k),'%04d'),'.tif'];
        temp_depth = imread(imname_depth);
        %         imname_depth = ['D1457_08Octopus_20220826_161130008_20221212_163306870_00000',...
        %             num2str(15512+imnum(k)),'_000000',num2str(imnum(k),'%04d'),...
        %             '_3D Depth Image (Orthographic to Reference).tiff'];
        %         temp_depth = imread([depthpath, filesep, imname_depth]);
        if k == 1
            depthmap = temp_depth(:,:,1);
        else
            depthmap(:,:,k) = temp_depth(:,:,1);
        end
    catch
    end
end

opt.cache = depthmap;
opt.cachenum = imnum;

% Flatten stack
depthmap = median(depthmap,3,'omitnan');

% Median filtering with special treatment of nans, to keep the size of
% valid data
% depthmap = nanmedfilt2(depthmap,[med2 med2]).*~isnan(depthmap);
temp = isnan(depthmap);
try
    hasgpu = ~isempty(gpuDeviceTable);
catch
    hasgpu = false;
end
if opt.fillmode=="fillmore"
    % Fill more, by setting larger median filter
    if hasgpu
        depthmap = nanmedfilt2gpu(depthmap,[31 31]);
    else
        depthmap = nanmedfilt2(depthmap,[31 31]);
    end
else
    % Use median filter settings
    if hasgpu
        depthmap = nanmedfilt2gpu(depthmap,[med2 med2]);
    else
        depthmap = nanmedfilt2(depthmap,[med2 med2]);
    end
end
if opt.fillmode=="nofill"
    depthmap(temp) = nan;
end

end