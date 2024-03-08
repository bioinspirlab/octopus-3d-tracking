%% Various analysis and plots for octopus data

function octo_ExportVideo_3Dviews(clip,conflictmode)
% Clip is the clip name
% Conflictmode is either "overwrite" or "skip", to do the same if the
% requested file already exists.

if nargin<1
clip = "O15_1611_7318_L1"; % Inspected with makima
clip = "O15_1611_7318_L2"; % Inspected with makima
clip = "O15_1611_7882_R1"; % Inspected with makima
clip = "O15_1611_7882_R3"; % Need rerun, added pt 37
clip = "O15_1611_7882_R4"; % Inspected with makima
clip = "O15_1611_11589_R3"; % Inspected with makima
clip = "O15_1611_11589_R4"; % Was rerun
% clip = "O15_1611_13209_L2"; % Inspected with makima
% clip = "O15_1611_13209_L3"; % Inspected with makima
clip = "O15_1611_15512_L1"; % Inspection looks good, quads generated, needs inspection
% clip = "O15_1611_19083_L1"; % Inspected with makima
clip = "O15_1611_19083_L2"; % Inspected with makima
% clip = "O15_1611_19083_L3"; % Need to rerun with contraints to spline
% clip = "O15_1611_22417_L1"; % Inspected with makima
% clip = "O15_1611_22417_L2"; % Needs inspecting
% clip = "O15_1611_22417_L3"; % Needs inspecting
% clip = "O15_1611_22417_L4"; % Needs inspecting
clip = "O14_24216_L3";
cliparr = {"O14_24216_L1","O14_24216_L2","O14_24216_L3","O10_7875_R1","O10_6946_L1","O15_1611_19083_L1",...
"O15_1611_19083_L2","O15_1611_19083_L3",...
"O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4",...
"O15_1611_7318_L1","O15_1611_7318_L2","O15_1611_7882_R1",...
"O15_1611_7882_R3","O15_1611_7882_R4","O15_1611_11589_R3",...
"O15_1611_11589_R4","O15_1611_13209_L2","O15_1611_13209_L3"};
clip = cliparr{16};
clip = "O15_1611_22417_L2";
end

if nargin<2
    conflictmode = "skip";
end


%% Data load

matpath = octo_InitializeData(clip);
load(matpath,'dropboxpath','vidpath');
preprocessdata = octo_PreProcess(clip);
load(preprocessdata,'xarr','yarr','arr','opt','ptdatmm','numframes','bgmove_mm')

%% Undo motion correction on points

ptdatmm(:,:,2:end) = ptdatmm(:,:,2:end)+bgmove_mm;
arr(:,:,2:end) = arr(:,:,2:end)+bgmove_mm;
% Spline, but not in mm but in pixels:
arrpix = arr; dx = median(diff(xarr)); dy = median(diff(yarr));
arrpix(:,1,:) = (arrpix(:,1,:)-xarr(1))/dx;
arrpix(:,2,:) = (arrpix(:,2,:)-yarr(1))/dy;


%% Make 3D plot video

outvid3d = [dropboxpath, filesep, 'visualizations', filesep ,...
    char(clip), '_', char(opt.method), '_3D.mp4'];

if exist(outvid3d,"file") && conflictmode~="overwrite"
    disp(['File ' outvid3d ' already exists. Skipping video export.']);
else

vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

vortho = VideoReader(vidpath);
[X,Y] = meshgrid(xarr,yarr);

ymin = min(ptdatmm(:,3,:),[],'all')-10;
ymax = max(ptdatmm(:,3,:),[],'all');
xmin = min(ptdatmm(:,1,:),[],'all')-10;
xmax = max(ptdatmm(:,1,:),[],'all')+10;
zmin = min(ptdatmm(:,2,:),[],'all')-10;
zmax = max(ptdatmm(:,2,:),[],'all')+10;
Z = X*0+ymin+10;

ff2 = figure(6);
ff2.Position = [100 100 1500 600];
for i=1:numframes
    fr = readFrame(vortho);

    % Initialize
    if i == 1
        % Set up axis 1
        subplot(1,2,1)
        bgimg1 = warp(X,Z,Y,fr);
        bgimg1.FaceAlpha = 0.7;
        hold on
        axis normal
        switch opt.method
            case "linear"
                pl3da = plot3(ptdatmm(:,1,i),ptdatmm(:,3,i),ptdatmm(:,2,i),...
                    'o-','LineWidth',3,'MarkerSize',5,'MarkerEdgeColor','r');
            case "spline"
                pl3da = plot3(ptdatmm(:,1,i),ptdatmm(:,3,i),ptdatmm(:,2,i),...
                    'o','MarkerSize',5,'MarkerEdgeColor','r');
                pl3darra = plot3(arr(:,1,i),arr(:,3,i),arr(:,2,i),...
                    '-b','LineWidth',3);
        end
        axis image
        axis manual
        ttta = title(['Frame: ', num2str(i,'%3.0f')]);
        pl3d1a = plot3(ptdatmm(1,1,i),ptdatmm(1,3,i),ptdatmm(1,2,i),...
            'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
        hold off
        xlabel('X (mm)')
        ylabel('Z (mm)')
        zlabel('Y (mm)')
        xlim([xmin,xmax]);
        zlim([zmin,zmax]);
        ylim([ymin,ymax])
        ax1 = gca;
        ax1.YDir = 'reverse';
        ax1.View = [17,12];      % Front view
        ax1.XGrid = "on";
        ax1.YGrid = "on";
        ax1.ZGrid = "on";

        % Set up axis 2
        subplot(1,2,2)
        bgimg2 = warp(X,Z,Y,fr);
        bgimg2.FaceAlpha = 0.7;
        hold on
        switch opt.method
            case "linear"
                pl3db = plot3(ptdatmm(:,1,i),ptdatmm(:,3,i),ptdatmm(:,2,i),...
                    'o-','LineWidth',3,'MarkerSize',5,'MarkerEdgeColor','r');
            case "spline"
                pl3db = plot3(ptdatmm(:,1,i),ptdatmm(:,3,i),ptdatmm(:,2,i),...
                    'o','MarkerSize',5,'MarkerEdgeColor','r');
                pl3darrb = plot3(arr(:,1,i),arr(:,3,i),arr(:,2,i),...
                    '-b','LineWidth',3);
        end
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
        ax2.XDir = 'reverse';
        ax2.View = [-100,-10];   % Side view
        ax2.XGrid = "on";
        ax2.YGrid = "on";
        ax2.ZGrid = "on";
    end

    % Update plots if this frame has any data
    if sum(~isnan(ptdatmm(:,:,i)),[1 2])>0

        bgimg1.CData = fr;
        bgimg2.CData = fr;
        switch opt.method
            case "linear"
            case "spline"
                pl3darra.XData = arr(:,1,i)';
                pl3darra.YData = arr(:,3,i)';
                pl3darra.ZData = arr(:,2,i)';
                pl3darrb.XData = arr(:,1,i)';
                pl3darrb.YData = arr(:,3,i)';
                pl3darrb.ZData = arr(:,2,i)';
        end
        pl3da.XData = ptdatmm(:,1,i)';
        pl3da.YData = ptdatmm(:,3,i)';
        pl3da.ZData = ptdatmm(:,2,i)';
        pl3db.XData = ptdatmm(:,1,i)';
        pl3db.YData = ptdatmm(:,3,i)';
        pl3db.ZData = ptdatmm(:,2,i)';
        pl3d1a.XData = ptdatmm(1,1,i)';
        pl3d1a.YData = ptdatmm(1,3,i)';
        pl3d1a.ZData = ptdatmm(1,2,i)';
        pl3d1b.XData = ptdatmm(1,1,i)';
        pl3d1b.YData = ptdatmm(1,3,i)';
        pl3d1b.ZData = ptdatmm(1,2,i)';
        ttta.String = ['Frame: ', num2str(i,'%3.0f')];
        tttb.String = ['Frame: ', num2str(i,'%3.0f')];

        drawnow;
        frame = getframe(ff2);
        frame = imresize(frame.cdata,2);
        writeVideo(vw,frame);
        %     pause(0.5)
    end
end
close(vw)
end
end
