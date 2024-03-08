%% Information

% Export simple depth maps as video and pngs. No points overlaid.

function octo_ExportVideoPNG_DepthMap(clip)
% clip is the clip name

%% Data load

matpath = octo_InitializeData(clip);
load(matpath,'dropboxpath','orthopath','startnum');
preprocessdata = octo_PreProcess(clip);
load(preprocessdata,'opt','numframes')

%% Make video of ortho depth map

outvid3d = [dropboxpath, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_depthmap.mp4'];
pngloc = [dropboxpath, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_png_depthmap', filesep];
[~,~] = mkdir(pngloc);
% pngloc = ['/Users/joost/Desktop/octo/clip_pngs/' char(clip) '_depthmap' filesep];
% pngloc = ['C:\Users\joost\Dropbox\MBARI\Projects\EyeRIS-octopus\images\' char(clip) '_depthmap' filesep];
vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

for i=1:numframes
    fpath = [orthopath num2str(startnum+i-1) '.tif'];
    fr = imread(fpath);
    fr = double(fr(:,:,1));

    if i==1
        f7 = figure(7);
        f7.Position(3) = 800;
        f7.Position(4) = 800;

        im1 = imagesc(fr);
        axis image
        ax1 = gca;
        cmap = colormap('turbo');
        % Reduce saturation
        cmapHSV = rgb2hsv(cmap);
        cmapHSV(:,2) = 0.7*cmapHSV(:,2);
        cmapHSV(:,3) = 0.95*cmapHSV(:,3);
        colormap(hsv2rgb(cmapHSV))
        clim([30 190])
%         clim([0 160])
    end
    im1.CData = fr;
drawnow;
% Grab overlaid frame
frame = getframe(ax1);
% frame = imresize(frame.cdata,[2000 2000]);
frame=frame.cdata;

% Save a PNG
imwrite(frame,[pngloc 'frame' num2str(i,'%04.0f') '.png'])

writeVideo(vw,frame);
%     pause(0.5)
%     if i>1510
%         pause(1)
%     end
end
close(vw)

