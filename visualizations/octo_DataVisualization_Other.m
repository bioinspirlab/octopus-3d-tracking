%% Various analysis and plots for octopus data

% Other visualization. Not quite up to date.

%% Inputs
clear all

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


%% Data load

% matpath = octo_InitializeData(clip);
% load(matpath);
% processdata = octo_PreProcess(clip);
% load(processdata)
matpath = octo_InitializeData(clip);
load(matpath);
corrdata = octo_DataCorrection(clip);
load(corrdata)
mode = opt.method;

stopherefornow  % Execute one of the section below for further analysis and visualization


%% Plot cumulative (euclidian) distance between points

% Calculate cumulative distance from point 1
newnumpoints = size(ptdatmm,1);
cumdist = nan(numframes,newnumpoints-1);

for i = 1:numframes
    temp = ~isnan(ptdatmm(:,1,i));
    temp = temp(:);
    for j = 2:newnumpoints
        temp1 = find(temp(1:j-1)>0,1,'last');
        if ~isempty(temp1)
            if j>2
                ref = cumdist(i,j-2);
            else
                ref = 0;
            end
            cumdist(i,j-1) = ref+sqrt((ptdatmm(j,1,i)-ptdatmm(temp1,1,i)).^2 + ...
                (ptdatmm(j,2,i)-ptdatmm(temp1,2,i)).^2 + ...
                (ptdatmm(j,3,i)-ptdatmm(temp1,3,i)).^2);
        end
    end
end


figure(4)
timax = (1:numframes)/framerate;
plot(timax,cumdist,'linewidth',2)
legend(legendCell = cellstr(num2str(ptorder(2:end)', '%-d')))
ylabel('Distance (mm)')
xlabel('Time (s)')
optimizeFig;


%% Make video of arm motion, with arm stretch and curvature in separate plots

mode = "stretch";
% mode = "curvature";

% Prepare data
outvid3d = [dropboxpath 'octo_' char(clip) '_3D_' char(mode) '.mp4'];
vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

[X,Y] = meshgrid(xarr,yarr);

% Get reference
refpos = ptdatmm(1,:,1);
ptdatmm = ptdatmm-refpos;
arr = arr-refpos;

margin = 20;
ymin = min(ptdatmm(:,3,:),[],'all')-margin;
ymax = max(ptdatmm(:,3,:),[],'all');
xmin = min(ptdatmm(:,1,:),[],'all')-margin;
xmax = max(ptdatmm(:,1,:),[],'all')+margin;
zmin = min(ptdatmm(:,2,:),[],'all')-margin;
zmax = max(ptdatmm(:,2,:),[],'all')+margin;
Z = X*0+ymin+10;


f9 = figure(9);
f9.Position = [173 799 1716 538];

% First draw arms
subplot(1,3,1);
i=1;
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
ax2.View = [-16,19];
ax2.XGrid = "on";
ax2.YGrid = "on";
ax2.ZGrid = "on";

% Load curve/stretch data
stretchfile = [char(clip) '_stretchcurve.mat'];
load(stretchfile,'curvax','cumdist','stretchax')

switch mode
    case "curvature"
        % Curvature
        subplot(1,3,2:3)
        plot(1:10)
        copyobj(allchild(curvax),gca);
        c = colorbar;
        c.Label.String = '(Bend radius)^{-1} (mm^{-1})';
        clim([0 20])
        cc=colormap;
        cc(1,:) = [0 0 0];
        colormap(cc)
        title('Arm curvature')

    case "stretch"
        % Stretch
        subplot(1,3,2:3)
        plot(1:10)
        copyobj(allchild(stretchax),gca);
        c = colorbar;
        c.Label.String = 'Stretch (% from median)';
        clim([-80 80])
        title('Arm stretch')
end


set(gca,'FontSize',15);
view([90 -90])
ylim([0,max(timax)])
xlim([0,1])
ylabel('Time (s)')
xlabel('Normalized arm position')


hold on
% Draw vertical line
wl = plot([0 1],[0 0],'-w','Linewidth',2);
hold off

% Loop through frames to visualize all data
for i=1:numframes
    pl3darrb.XData = arr(:,1,i)';
    pl3darrb.YData = arr(:,3,i)';
    pl3darrb.ZData = arr(:,2,i)';
    if numel(guidepts)>0
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

    pl3d1b.XData = ptdatmm(1,1,i)';
    pl3d1b.YData = ptdatmm(1,3,i)';
    pl3d1b.ZData = ptdatmm(1,2,i)';
    tttb.String = ['Frame: ', num2str(i,'%3.0f')];

    wl.YData = [timax(i) timax(i)];

    drawnow;
    % Grab 3D spline plot
    frame = getframe(f9);
%     frame = imresize(frame.cdata,(frH/2)/size(frame.cdata,1));
%     maxW = min(size(frame,2),frW-frH);
%     outfr(1:frH/2,frH+1:frH+maxW,:) = frame(1:frH/2,1:maxW,:);
    writeVideo(vw,frame);
end
close(vw)


clear ptdatmm % Because we have messed it up
clear arr

%% %% Make video of arm motion, colored by curvature


% Prepare data
outvid3d = [dropboxpath 'octo_' char(clip) '_3D_colorcurvature.mp4'];
vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

[X,Y] = meshgrid(xarr,yarr);

% Load curve/stretch data
stretchfile = [char(clip) '_stretchcurve.mat'];
load(stretchfile,'curvax','cumdist','splpts','curvdat')

% Get reference
refpos = ptdatmm(1,:,1);
ptdatmm = ptdatmm-refpos;
arr = arr-refpos;
splpts = splpts - refpos;

margin = 20;
ymin = min(ptdatmm(:,3,:),[],'all')-margin;
ymax = max(ptdatmm(:,3,:),[],'all');
xmin = min(ptdatmm(:,1,:),[],'all')-margin;
xmax = max(ptdatmm(:,1,:),[],'all')+margin;
zmin = min(ptdatmm(:,2,:),[],'all')-margin;
zmax = max(ptdatmm(:,2,:),[],'all')+margin;
Z = X*0+ymin+10;



f9 = figure(9);
% f9.Position = [173 799 1716 538];

minval = 0;
maxval = 15;
clrscale = jet((maxval-minval)*10+1);

% First draw arms
i=1;
pl3db = plot3(ptdatmm(~guideptind,1,i),ptdatmm(~guideptind,3,i),ptdatmm(~guideptind,2,i),...
    'o','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on
if numel(guidepts)>0
    pl3dbGP = plot3(ptdatmm(guideptind,1,i),ptdatmm(guideptind,3,i),ptdatmm(guideptind,2,i),...
        'o','MarkerSize',7,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);
end
for j=1:size(splpts,1)-1
    val = min(max(curvdat(i,j),minval),maxval);
    clrt = clrscale(round(val*10)+1,:);
pl3darrb(j) = plot3(splpts(j,1,i),splpts(j,3,i),splpts(j,2,i),...
    '.','Color',clrt,'MarkerSize',10);
end
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
% ax2.View = [17,12];      % Front view
ax2.View = [-12,40];
ax2.XGrid = "on";
ax2.YGrid = "on";
ax2.ZGrid = "on";

cb = colorbar;
caxis([minval maxval])
cb.Label.String = '(Bend radius)^{-1} (mm^{-1})';
ax2.FontSize = 15;

% Loop through frames to visualize all data
for i=1:numframes
%     pl3darrb.XData = arr(:,1,i)';
%     pl3darrb.YData = arr(:,3,i)';
%     pl3darrb.ZData = arr(:,2,i)';
    for j=1:size(splpts,1)-1
        val = min(max(curvdat(i,j),minval),maxval);
        clrt = clrscale(round(val*10)+1,:);
%         pl3darrb(j) = plot3(splpts(:,1,i),splpts(:,3,i),splpts(:,2,i),...
%             '.','Color',clrt);
        pl3darrb(j).XData = splpts(j,1,i);
        pl3darrb(j).YData = splpts(j,3,i);
        pl3darrb(j).ZData = splpts(j,2,i);
        pl3darrb(j).Color = clrt;
    end

    if numel(guidepts)>0
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

    pl3d1b.XData = ptdatmm(1,1,i)';
    pl3d1b.YData = ptdatmm(1,3,i)';
    pl3d1b.ZData = ptdatmm(1,2,i)';
    tttb.String = ['Frame: ', num2str(i,'%3.0f')];

    drawnow;
    % Grab 3D spline plot
    frame = getframe(f9);
%     frame = imresize(frame.cdata,(frH/2)/size(frame.cdata,1));
%     maxW = min(size(frame,2),frW-frH);
%     outfr(1:frH/2,frH+1:frH+maxW,:) = frame(1:frH/2,1:maxW,:);
    writeVideo(vw,frame);
end
close(vw)


clear ptdatmm % Because we have messed it up
clear arr


%% Plot multiple arms at different times before and after release
clear all

prereleaseframe = -120;
postreleaseframe = 120;
precolor = [0 0.447 0.741];
releasecolor = [0.85 0.325 0.098];
postcolor = [0.929 0.694 0.125];
binint = 0.05;
ymax = 14;
ymin = 0;

figure(11)

% L2
subplot(4,2,3)
clip = "O15_1611_22417_L2";
attachloc = [0.5398, 0.7027];
touchframe = 950;
load(['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'],'timax','curvdat','cumdist','yyy','ptdistarr');
xarr = cumdist(1:end-1)/max(cumdist);
[N,~,bin] = histcounts(xarr,0:binint:1);
numbins = numel(N);
xarr2 = [];
curvpre = [];
curvduring = [];
curvpost = [];
for i=1:numbins
    xarr2(i) = mean(xarr(bin==i));
    curvpre(i) = mean(curvdat(touchframe+prereleaseframe,bin==i));
    curvduring(i) = mean(curvdat(touchframe,bin==i));
    curvpost(i) = mean(curvdat(touchframe+postreleaseframe,bin==i));
end
plot(xarr2,curvpre,'Color',precolor)
hold on
plot(xarr2,curvduring,'Color',releasecolor,'Linewidth',3)
plot(xarr2,curvpost,'Color',postcolor)
% Draw touch area
ylim([ymin ymax]);
yl = ylim;
rectangle('Position',[attachloc(1) yl(1)+0.01*diff(yl) diff(attachloc) diff(yl)*0.05],...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
hold off
hax = gca;
hax.Children = circshift(hax.Children, -1);
ylabel('(Bend radius)^{-1} (mm^{-1})')
title('L2')
legend('Prerelease','Release','Postrelease')


% L3
subplot(4,2,5)
clip = "O15_1611_19083_L3";
attachloc = [0.58, 0.72];
touchframe = 2911;
load(['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'],'timax','curvdat','cumdist');
xarr = cumdist(1:end-1)/max(cumdist);
[N,~,bin] = histcounts(xarr,0:binint:1);
numbins = numel(N);
xarr2 = [];
curvpre = [];
curvduring = [];
curvpost = [];
for i=1:numbins
    xarr2(i) = mean(xarr(bin==i));
    curvpre(i) = mean(curvdat(touchframe+prereleaseframe,bin==i));
    curvduring(i) = mean(curvdat(touchframe,bin==i));
    curvpost(i) = mean(curvdat(touchframe+postreleaseframe,bin==i));
end
plot(xarr2,curvpre,'Color',precolor)
hold on
plot(xarr2,curvduring,'Color',releasecolor,'Linewidth',3)
plot(xarr2,curvpost,'Color',postcolor)
% Draw touch area
ylim([ymin ymax]);
yl = ylim;
rectangle('Position',[attachloc(1) yl(1)+0.01*diff(yl) diff(attachloc) diff(yl)*0.05],...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
hold off
hax = gca;
hax.Children = circshift(hax.Children, -1);
ylabel('(Bend radius)^{-1} (mm^{-1})')
title('L3')
legend('Prerelease','Release','Postrelease');

% R3
subplot(4,2,6)
clip = "O15_1611_7882_R3";
% attachloc = [0.3, 0.4];
touchframe = 1561;
load(['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'],'timax','curvdat','cumdist');
xarr = cumdist(1:end-1)/max(cumdist);
[N,~,bin] = histcounts(xarr,0:binint:1);
numbins = numel(N);
xarr2 = [];
curvpre = [];
curvduring = [];
curvpost = [];
for i=1:numbins
    xarr2(i) = mean(xarr(bin==i));
    curvpre(i) = mean(curvdat(touchframe+prereleaseframe,bin==i));
    curvduring(i) = mean(curvdat(touchframe,bin==i));
    curvpost(i) = mean(curvdat(touchframe+postreleaseframe,bin==i));
end
plot(xarr2,curvpre,'Color',precolor)
hold on
plot(xarr2,curvduring,'Color',releasecolor,'Linewidth',3)
plot(xarr2,curvpost,'Color',postcolor)
% Draw touch area
ylim([ymin ymax]);
yl = ylim;
rectangle('Position',[attachloc(1) yl(1)+0.01*diff(yl) diff(attachloc) diff(yl)*0.05],...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
hold off
hax = gca;
hax.Children = circshift(hax.Children, -1);
ylabel('(Bend radius)^{-1} (mm^{-1})')
title('R3')
legend('Prerelease','Release','Postrelease')

% L4
subplot(4,2,7)
clip = "O15_1611_22417_L4";
touchframe = 1697;
load(['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'],'timax','curvdat','cumdist');
xarr = cumdist(1:end-1)/max(cumdist);
[N,~,bin] = histcounts(xarr,0:binint:1);
numbins = numel(N);
xarr2 = [];
curvpre = [];
curvduring = [];
curvpost = [];
for i=1:numbins
    xarr2(i) = mean(xarr(bin==i));
    curvpre(i) = mean(curvdat(touchframe+prereleaseframe,bin==i));
    curvduring(i) = mean(curvdat(touchframe,bin==i));
    curvpost(i) = mean(curvdat(touchframe+postreleaseframe,bin==i));
end
plot(xarr2,curvpre,'Color',precolor)
hold on
plot(xarr2,curvduring,'Color',releasecolor,'Linewidth',3)
plot(xarr2,curvpost,'Color',postcolor)
% Draw touch area
ylim([ymin ymax]);
yl = ylim;
rectangle('Position',[attachloc(1) yl(1)+0.01*diff(yl) diff(attachloc) diff(yl)*0.05],...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
hold off
hax = gca;
hax.Children = circshift(hax.Children, -1);
ylabel('(Bend radius)^{-1} (mm^{-1})')
title('L4')
legend('Prerelease','Release','Postrelease')


% R4
subplot(4,2,7)
clip = "O15_1611_11589_R4";
touchframe = 906;
load(['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'],'timax','curvdat','cumdist');
xarr = cumdist(1:end-1)/max(cumdist);
[N,~,bin] = histcounts(xarr,0:binint:1);
numbins = numel(N);
xarr2 = [];
curvpre = [];
curvduring = [];
curvpost = [];
for i=1:numbins
    xarr2(i) = mean(xarr(bin==i));
    curvpre(i) = mean(curvdat(touchframe+prereleaseframe,bin==i));
    curvduring(i) = mean(curvdat(touchframe,bin==i));
    curvpost(i) = mean(curvdat(touchframe+postreleaseframe,bin==i));
end
plot(xarr2,curvpre,'Color',precolor)
hold on
plot(xarr2,curvduring,'Color',releasecolor,'Linewidth',3)
plot(xarr2,curvpost,'Color',postcolor)
% Draw touch area
ylim([ymin ymax]);
yl = ylim;
rectangle('Position',[attachloc(1) yl(1)+0.01*diff(yl) diff(attachloc) diff(yl)*0.05],...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
hold off
hax = gca;
hax.Children = circshift(hax.Children, -1);
ylabel('(Bend radius)^{-1} (mm^{-1})')
title('R4')
legend('Prerelease','Release','Postrelease')

%% Make 3D view from PLY
% The PLY needs to be exported with texture from RXlive, and then opened
% and saved in MeshLab, which corrects invalid values.

plypath = '/Users/joost/Dropbox/MBARI/Projects/EyeRIS-octopus/D1457_08Octopus_20220826_161130008_20230928_231216883_0000017877_0000000000_3D_Mesh_mod.ply';

% viewer = viewer3d;
% viewer.BackgroundColor = "white";

meshvar = readSurfaceMesh(plypath);
sshow = surfaceMeshShow(meshvar,BackgroundColor="white");
