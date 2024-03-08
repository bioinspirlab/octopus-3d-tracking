function outvid = octo_ExportVideoOverlayPointsAndSpline(clip,conflictmode)
% Clip is the clip name
% Conflictmode is either "overwrite" or "skip", to do the same if the
% requested file already exists.


%% Load Data
matpath = octo_InitializeData(clip);
load(matpath,'dropboxpath','framerate','vidpath','orthopath','startnum');
preprocessdata = octo_PreProcess(clip);
load(preprocessdata,'ptdat','ptdatmm','arr','xarr','yarr','bgmove_px',...
    'bgmove_mm','ptnames','guidepts','opt')

%% Set settings


%% Initialize

outvid = [dropboxpath, filesep, 'visualizations', filesep,...
    char(clip), '_', opt.method, '_pointoverlay.mp4'];
if nargin<2
    conflictmode = "skip";
end

if exist(outvid,"file") && conflictmode~="overwrite"
    disp(['File ' outvid ' already exists. Skipping video export.']);
else

% Open ortho video
vortho = VideoReader(vidpath);
if ~exist('framerate','var')
    framerate = vortho.FrameRate;
end

% Initialize video
vw = VideoWriter(outvid,'MPEG-4');
vw.Quality = 95;
open(vw);

numframes = size(ptdat,1);
numpoints = size(ptdat,3);

%% Undo motion correction on points
ptdat(2:end,:,:) = ptdat(2:end,:,:)+bgmove_px;
ptdatmm(:,:,2:end) = ptdatmm(:,:,2:end)+bgmove_mm;
arr(:,:,2:end) = arr(:,:,2:end)+bgmove_mm;
% Remove spline points where no z coordinate is available
arr(:,1,:) = arr(:,1,:).*~isnan(arr(:,3,:));
arr(:,2,:) = arr(:,2,:).*~isnan(arr(:,3,:));
% Spline, but not in mm but in pixels:
arrpix = arr; dx = median(diff(xarr)); dy = median(diff(yarr));
arrpix(:,1,:) = (arrpix(:,1,:)-xarr(1))/dx;
arrpix(:,2,:) = (arrpix(:,2,:)-yarr(1))/dy;

%% Loop through data and frames
firstloop = true;
for i=1:numframes
    temp_ortho = readFrame(vortho);

    % Only continue if there is data for this frame
    if sum(~isnan(sum(ptdat(i,:,:),2)))>0

        imname_depth = [orthopath,num2str(startnum-1+i,'%04d'),'.tif'];
        temp_depth = imread(imname_depth);

        % Initialize figure position
        if firstloop
            ff = figure(3);
            ff.Name = 'Octo_ortho';
            ff.Position(2) = 100;
            ff.Position(3) = 1.2*size(temp_ortho,2);
            ff.Position(4) = 0.6*size(temp_ortho,1);
            scrsz = get(0,'ScreenSize');
            if scrsz(3) < ff.Position(1)+ff.Position(3)
                ff.Position(1) = scrsz(3)-ff.Position(3)-60;
            end
            if scrsz(4) < ff.Position(2)+ff.Position(4)
                ff.Position(2) = scrsz(4)-ff.Position(4)-60;
            end
            
        end

        [depthmap,opt] = octo_ProcessDepthMap(i,orthopath,startnum,opt);

        xloc = ptdat(i,1,:);
        yloc = ptdat(i,2,:);
%         ptdatmm(:,1,i) = permute(interp1(xpix,xarr,xloc),[3,1,2]);
%         ptdatmm(:,2,i) = permute(interp1(ypix,yarr,yloc),[3,1,2]);
%         ptdatmm(:,3,i) = permute(interp2(X,Y,depthmap,xloc,yloc),[3,1,2]);

        if firstloop
            figure(ff)
            subplot(1,2,1)
            immL = imagesc(temp_ortho(:,:,1));
            colormap('gray')
            hold on
            pp1 = plot(reshape(ptdat(i,1,:),[numpoints 1]),...
                reshape(ptdat(i,2,:),[numpoints 1]),'ro','MarkerSize',8);
            ak1 = plot(1:5,1:5);
            % Add numbers of points
            for k = 1:numpoints
                if sum(k==guidepts)>0
                    % This is a guide point (white color)
                    clr = 'white';
                else
                    clr = 'red';
                end
                tx1(k) = text(ptdat(i,1,k),ptdat(i,2,k),['  ' ptnames{k}],...
                    'HorizontalAlignment','left','Color',clr,'FontSize',11);
            end
            hold off
            ax1 = gca;
            axis image
            ax1.YDir = 'reverse';
            figtxt1 = text(ax1.XLim(1)+10,ax1.YLim(2)-10, 'Frame: 1',...
                'HorizontalAlignment','left','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',12,'Color','white');

            subplot(1,2,2)
            immR = imagesc(xarr,yarr,depthmap);
            ax2 = gca;
            colormap(ax2,"turbo");
            hold on
%             pp2 = plot(reshape(ptdat(i,1,:),[numpoints 1]),...
%                 reshape(ptdat(i,2,:),[numpoints 1]),'ro','MarkerSize',8);
            pp2 = plot(reshape(ptdatmm(:,1,i),[numpoints 1]),...
                reshape(ptdatmm(:,2,i),[numpoints 1]),'ro','MarkerSize',8);
            ak = plot(1:5,1:5);
            % Add numbers
            for k = 1:numpoints
                if sum(k==guidepts)>0
                    % This is a guide point (white color)
                    clr = 'white';
                else
                    clr = 'red';
                end
                tx2(k) = text(ptdatmm(k,1,i),ptdatmm(k,2,i),['  ' ptnames{k}],...
                    'HorizontalAlignment','left','Color',clr,'FontSize',11);
            end
            hold off

            axis image
%             ax2.YDir = 'reverse';
            ax2.YDir = "normal";
            figtxt2 = text(xarr(10),yarr(end-10), 'Frame: 1',...
                'HorizontalAlignment','left','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',12,'Color','white');
            firstloop = false;
        end

        % Update image and spline data
        immL.CData = temp_ortho(:,:,1);
        immR.CData = depthmap;

        % Update point data
        pp1.XData = reshape(ptdat(i,1,:),[numpoints 1]);
        pp1.YData = reshape(ptdat(i,2,:),[numpoints 1]);
        pp2.XData = reshape(ptdatmm(:,1,i),[numpoints 1]);
        pp2.YData = reshape(ptdatmm(:,2,i),[numpoints 1]);

        % Update label positions
        for j = 1:numpoints
            tx1(j).Position = [ptdat(i,1,j), ptdat(i,2,j)];
            tx2(j).Position = [ptdatmm(j,1,i), ptdatmm(j,2,i)];
        end

        % Update spline data
        ak1.XData = arrpix(:,1,i);
        ak1.YData = arrpix(:,2,i);
        ak.XData = arr(:,1,i);
        ak.YData = arr(:,2,i);

        % Update frame label
        figtxt1.String = ['Frame: ', num2str(i,'%3.0f')];
        figtxt2.String = ['Frame: ', num2str(i,'%3.0f')];

        pause(0.05);
        try
            frame1 = getframe(ax1);
            frame2 = getframe(ax2);
            frame = zeros(max([size(frame1.cdata);size(frame2.cdata)]).*[1,2,1],'uint8');
            frame(1:size(frame1.cdata,1),1:size(frame1.cdata,2),:) = frame1.cdata;
            xref = size(frame,2)/2;
            frame((1:size(frame2.cdata,1)),xref+(1:size(frame2.cdata,2)),:) = frame2.cdata;
            frame = imresize(frame,2);
            writeVideo(vw,frame);
        catch
            warning(['Frame ' num2str(i) ' not captured to video']);
        end
    end
end
close(vw)
end
end