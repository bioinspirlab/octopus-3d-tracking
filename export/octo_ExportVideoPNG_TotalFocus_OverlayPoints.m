%% Information

% Export processed images for given clip(s). This is different from
% octo_ExportVideoOverlayPointsAndSpline which does not allow multiple
% clips to be spliced, and gives a two-panel output.

function octo_ExportVideoPNG_TotalFocus_OverlayPoints(indata)
% indata is the clip name; if multiple are given, they should be different
% arm datasets on the same clip

%% Data load

matpath = octo_InitializeData(indata(1));
load(matpath,'dropboxpath','vidpath');
preprocessdata = octo_PreProcess(indata(1));
load(preprocessdata,'opt','numframes','xarr','yarr')

%% Plot just grayscale total focus image with points and spline overlay

% multiarm = ["O15_1611_22417_L1","O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4"];
% multiarm = ["O15_1611_11589_R3","O15_1611_11589_R4"];
% multiarm = ["O15_1611_15512_L1"];
% indata = ["O15_1611_19083_L1","O15_1611_19083_L2","O15_1611_19083_L3"];
multiarm = indata;
for j=1:numel(multiarm)
    clip=multiarm(j);
    preprocessdata = octo_PreProcess(clip);
    S.(multiarm(j)) = load(preprocessdata,'ptdatmm','arr','xarr','yarr',...
        'guidepts','guideptind','bgmove_mm','numframes');

    % Undo background motion removal
S.(multiarm(j)).ptdatmm(:,:,2:end) = S.(multiarm(j)).ptdatmm(:,:,2:end)+S.(multiarm(j)).bgmove_mm;
S.(multiarm(j)).arr(:,:,2:end) = S.(multiarm(j)).arr(:,:,2:end)+S.(multiarm(j)).bgmove_mm;
end

outvid3d = [dropboxpath, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_totalfocus_ptoverlay.mp4'];
pngloc = [dropboxpath, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_png_totalfocus', filesep];
[status,msg] = mkdir(pngloc);
% pngloc = ['/Users/joost/Desktop/octo/clip_pngs/' char(clip) filesep];
vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

vortho = VideoReader(vidpath);
% [X,Y] = meshgrid(xarr,yarr);

% vfocusdepth = VideoReader(focusdepthpath);

% vquad = VideoReader(quadpath);
% vquad.CurrentTime = quadstartframe/vquad.FrameRate;
% quadfrnum = 0;

% Frame info
frH = 1080;
frW = 1920;

outfr = zeros(frH,frW,3,'uint8');
image(outfr)
axis image

markertype = ["-b","--b",":b","-.b"];
if numel(multiarm)>4
    markertype(5:numel(multiarm)) = "-b";
end

for i=1:numframes

    fr = readFrame(vortho);
    if i == 1
        f7 = figure(7);
        f7.Position(3) = 800;
        f7.Position(4) = 800;

        bgimg1 = image(xarr,yarr,fr);
        hold on

for j = 1:numel(multiarm)
        pl3da(j) = plot(S.(multiarm(j)).ptdatmm(~S.(multiarm(j)).guideptind,1,i),S.(multiarm(j)).ptdatmm(~S.(multiarm(j)).guideptind,2,i),...
            'o','MarkerSize',18,'MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        if numel(S.(multiarm(j)).guidepts)>0
            pl3daGP(j) = plot(S.(multiarm(j)).ptdatmm(S.(multiarm(j)).guideptind,1,i),S.(multiarm(j)).ptdatmm(S.(multiarm(j)).guideptind,2,i),...
                'o','MarkerSize',18,'MarkerEdgeColor','r','MarkerFaceColor','none');
        end
        pl3darra(j) = plot(S.(multiarm(j)).arr(:,1,i),S.(multiarm(j)).arr(:,2,i),...
            markertype(j),'LineWidth',7);
        hold on
        
        axis image
        axis manual
        %         ttta = title(['Frame: ', num2str(i,'%3.0f')]);
        
%         pl3d1a(j) = plot(S.(multiarm(j)).ptdatmm(1,1,i),S.(multiarm(j)).ptdatmm(1,2,i),...
%             'o','MarkerSize',26,'MarkerEdgeColor','r','MarkerFaceColor','r');
end
ttta = text(min(xarr)+5,min(yarr)+5,['t = ' num2str(i/vortho.FrameRate,'%03.2f') ' s'],...
            'VerticalAlignment','bottom','HorizontalAlignment','left',...
            'Color','white','FontSize',26,'Fontweight','bold','BackgroundColor',[0.5 0.5 0.5]);
        hold off
        xlabel('X (mm)')
        ylabel('Y (mm)')
        ax1 = gca;
        ax1.YDir = 'normal';
    end

% bgimg1.CData = fr;
% bgimg1.CData = uint8(3*double(fr).^1.1);
% bgimg1.CData = uint8(3*double(fr).^0.9);% USe this one for clip 11589, 15512
bgimg1.CData = uint8(4*double(fr).^0.9);% USe this one for clip 11589, 15512
% bgimg1.CData = uint8(6*double(fr).^0.9);% USe this one for clip 22417

% bgimg1.CData = 4.5*fr;    % The multiplier makes the image brighter
%     bgimg2.CData = fr;
for j = 1:numel(multiarm)
        pl3darra(j).XData = S.(multiarm(j)).arr(:,1,i)';
        pl3darra(j).YData = S.(multiarm(j)).arr(:,2,i)';
pl3da(j).XData = S.(multiarm(j)).ptdatmm(~S.(multiarm(j)).guideptind,1,i)';
pl3da(j).YData = S.(multiarm(j)).ptdatmm(~S.(multiarm(j)).guideptind,2,i)';
if numel(S.(multiarm(j)).guidepts)>0
pl3daGP(j).XData = S.(multiarm(j)).ptdatmm(S.(multiarm(j)).guideptind,1,i)';
pl3daGP(j).YData = S.(multiarm(j)).ptdatmm(S.(multiarm(j)).guideptind,2,i)';
end


% pl3d1a(j).XData = S.(multiarm(j)).ptdatmm(1,1,i)';
% pl3d1a(j).YData = S.(multiarm(j)).ptdatmm(1,2,i)';
end
% ttta.String = ['t = ' num2str(i/vortho.FrameRate,'%03.2f') ' s'];
ttta.String = [''];

drawnow;
% Grab overlaid frame
frame = getframe(ax1);
% frame = imresize(frame.cdata,[2000 2000]);
frame=frame.cdata;

% Save a PNG
imwrite(frame,[pngloc 'frame' num2str(i,'%04.0f') '.png'])

writeVideo(vw,frame);
end
close(vw)
end