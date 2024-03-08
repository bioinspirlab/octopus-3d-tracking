function octo_ExportVideo_Quad3D(clip,conflictmode)
% Clip is the clip name
% Conflictmode is either "overwrite" or "skip", to do the same if the
% requested file already exists.

if nargin<1
    cliparr = {"O14_24216_L1","O14_24216_L2","O14_24216_L3","O10_7875_R1","O10_6946_L1","O15_1611_19083_L1",...
        "O15_1611_19083_L2","O15_1611_19083_L3",...
        "O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4",...
        "O15_1611_7318_L1","O15_1611_7318_L2","O15_1611_7882_R1",...
        "O15_1611_7882_R3","O15_1611_7882_R4","O15_1611_11589_R3",...
        "O15_1611_11589_R4","O15_1611_13209_L2","O15_1611_13209_L3"};
    clip = cliparr{2};
end
if nargin<2
    conflictmode = "skip";
end


%% Data load

matpath = octo_InitializeData(clip);
load(matpath,'dropboxpath','vidpath','focusdepthpath','quadpath','quadstartframe');
preprocessdata = octo_PreProcess(clip);
load(preprocessdata,'xarr','yarr','arr','opt','ptdatmm','numframes',...
    'bgmove_mm','guideptind','guidepts','ptorder','ptnames')


%% Undo motion correction on points

ptdatmm(:,:,2:end) = ptdatmm(:,:,2:end)+bgmove_mm;
arr(:,:,2:end) = arr(:,:,2:end)+bgmove_mm;
% Spline, but not in mm but in pixels:
arrpix = arr; dx = median(diff(xarr)); dy = median(diff(yarr));
arrpix(:,1,:) = (arrpix(:,1,:)-xarr(1))/dx;
arrpix(:,2,:) = (arrpix(:,2,:)-yarr(1))/dy;


%% Make full quad

outvid3d = [dropboxpath, filesep, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_quad.mp4'];

if exist(outvid3d,"file") && conflictmode~="overwrite"
    disp(['File ' outvid3d ' already exists. Skipping video export.']);
else

vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

vortho = VideoReader(vidpath);
[X,Y] = meshgrid(xarr,yarr);

vfocusdepth = VideoReader(focusdepthpath);

vquad = VideoReader(quadpath);
vquad.CurrentTime = quadstartframe/vquad.FrameRate;
quadfrnum = 0;

% Frame info
frH = 1080;
frW = 1920;

outfr = zeros(frH,frW,3,'uint8');
image(outfr)
axis image

margin = 10;
ymin = min(arr(:,3,:),[],'all')-margin;
ymax = max(arr(:,3,:),[],'all');
xmin = min(arr(:,1,:),[],'all')-margin;
xmax = max(arr(:,1,:),[],'all')+margin;
zmin = min(arr(:,2,:),[],'all')-margin;
zmax = max(arr(:,2,:),[],'all')+margin;
Z = X*0+ymin+10;
% 
% ff2 = figure(7);
% ff2.Position(3) = 1500;
% ff2.Position(4) = 600;


for i=1:numframes

% Read quad frame and use bottom half
if i*vquad.FrameRate/60>quadfrnum
    quadfr = readFrame(vquad);
    outfr(frH/2:end,:,:) = quadfr(frH/2:end,:,:);
    quadfrnum = quadfrnum+1;
end

% Read focused depth map and put in top left
focdepthfr = readFrame(vfocusdepth);
outfr(1:frH/2,1:frH/2,:) = imresize(focdepthfr,[frH/2,frH/2]);


fr = readFrame(vortho);
if i == 1
    f7 = figure(7);
    f7.Position(3) = (frH/2);
    f7.Position(4) = frH/2;

    bgimg1 = image(xarr,yarr,fr);
    hold on
            pl3da = plot(ptdatmm(~guideptind,1,i),ptdatmm(~guideptind,2,i),...
                'o','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');
            hold on
            if numel(guidepts)>0
            pl3daGP = plot(ptdatmm(guideptind,1,i),ptdatmm(guideptind,2,i),...
                'o','MarkerSize',8,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);
            end
            pl3darra = plot(arr(:,1,i),arr(:,2,i),...
                '-b','LineWidth',3);
    hold on
    for jk = 1:numel(ptorder)
        pttxt(jk) = text(ptdatmm(jk,1,i),ptdatmm(jk,2,i),ptnames(jk),...
            'Color','r','VerticalAlignment','baseline','FontSize',12);
    end
    axis image
    axis manual
    %         ttta = title(['Frame: ', num2str(i,'%3.0f')]);
    ttta = text(min(xarr),min(yarr),['Frame: ', num2str(i,'%3.0f')],...
        'VerticalAlignment','bottom','HorizontalAlignment','left',...
        'Color','white','FontSize',12);
    pl3d1a = plot(ptdatmm(1,1,i),ptdatmm(1,2,i),...
        'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');

    hold off
    xlabel('X (mm)')
    ylabel('Y (mm)')
    ax1 = gca;
    ax1.YDir = 'normal';


    % Set up axis 2
    f8 = figure(8);
    f8.Position(3) = (frW-frH);
    f8.Position(4) = frH/2;
%             pl3db = plot3(ptdatmm(:,1,i),ptdatmm(:,3,i),ptdatmm(:,2,i),...
%                 'o','MarkerSize',7,'MarkerEdgeColor','r');
            pl3db = plot3(ptdatmm(~guideptind,1,i),ptdatmm(~guideptind,3,i),ptdatmm(~guideptind,2,i),...
                'o','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','r');
            hold on
            if numel(guidepts)>0
            pl3dbGP = plot3(ptdatmm(guideptind,1,i),ptdatmm(guideptind,3,i),ptdatmm(guideptind,2,i),...
                'o','MarkerSize',7,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);
            end
            pl3darrb = plot3(arr(:,1,i),arr(:,3,i),arr(:,2,i),...
                '-b','LineWidth',3);
    hold on
    axis image
    axis manual
    tttb = title(['Frame: ', num2str(i,'%3.0f')]);
    pl3d1b = plot3(ptdatmm(1,1,i),ptdatmm(1,3,i),ptdatmm(1,2,i),...
        'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
    hold off
    xlabel('X (mm)')
    ylabel('Z (mm)')
    zlabel('Y (mm)')
    xlim([xmin, xmax]);
    zlim([zmin,zmax]);
    ylim([ymin,ymax])
    ax2 = gca;
    ax2.YDir = 'reverse';
    ax2.View = [17,12];      % Front view
%     ax2.View = [-87,44];      % Other view
    ax2.XGrid = "on";
    ax2.YGrid = "on";
    ax2.ZGrid = "on";
end

bgimg1.CData = 30+fr;       % Just a little brighter
% bgimg1.CData = fr;    % Normal
%     bgimg2.CData = fr;

        pl3darra.XData = arr(:,1,i)';
        pl3darra.YData = arr(:,2,i)';
        %             pl3darra.ZData = arr(:,2,i)';
        pl3darrb.XData = arr(:,1,i)';
        pl3darrb.YData = arr(:,3,i)';
        pl3darrb.ZData = arr(:,2,i)';
pl3da.XData = ptdatmm(~guideptind,1,i)';
pl3da.YData = ptdatmm(~guideptind,2,i)';
if numel(guidepts)>0
pl3daGP.XData = ptdatmm(guideptind,1,i)';
pl3daGP.YData = ptdatmm(guideptind,2,i)';
pl3dbGP.XData = ptdatmm(guideptind,1,i)';
pl3dbGP.YData = ptdatmm(guideptind,3,i)';
pl3dbGP.ZData = ptdatmm(guideptind,2,i)';
end
pl3db.XData = ptdatmm(~guideptind,1,i)';
pl3db.YData = ptdatmm(~guideptind,3,i)';
pl3db.ZData = ptdatmm(~guideptind,2,i)';


for jk = 1:numel(ptorder)
    pttxt(jk).Position = [ptdatmm(jk,1,i),ptdatmm(jk,2,i)];
end

pl3d1a.XData = ptdatmm(1,1,i)';
pl3d1a.YData = ptdatmm(1,2,i)';
%     pl3d1a.ZData = ptdatmm(1,2,i)';
pl3d1b.XData = ptdatmm(1,1,i)';
pl3d1b.YData = ptdatmm(1,3,i)';
pl3d1b.ZData = ptdatmm(1,2,i)';
ttta.String = ['Frame: ', num2str(i,'%3.0f')];
tttb.String = ['Frame: ', num2str(i,'%3.0f')];

drawnow;
% Grab overlaid frame
frame = getframe(ax1);
frame = imresize(frame.cdata,(frH/2)/size(frame.cdata,1));
maxW = min(size(frame,2),frH/2);
outfr(1:frH/2,frH/2+1:frH/2+maxW,:) = frame(1:frH/2,1:maxW,:);
% Grab 3D spline plot
frame = getframe(f8);
frame = imresize(frame.cdata,(frH/2)/size(frame.cdata,1));
maxW = min(size(frame,2),frW-frH);
outfr(1:frH/2,frH+1:frH+maxW,:) = frame(1:frH/2,1:maxW,:);
writeVideo(vw,outfr);
%     pause(0.5)
%     if i>1510
%         pause(1)
%     end
end
close(vw)
end
end