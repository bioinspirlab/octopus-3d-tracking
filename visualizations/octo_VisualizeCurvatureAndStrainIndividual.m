%% Various analysis and plots for octopus data

if exist("clip")~=1
    error("Define a clip to visualize.")
end
if exist("input_data_path")~=1
    error("Define a data location path.")
end
if exist("output_data_path")~=1
    error("Define an output path.")
end

%% Settings

switch clip
    case "O10_6946_L1"   % Verified. Note: use adjusted reference, option (4) below
    case "O10_7875_R1"  % Verified. Clip too short for stretch reference
    case "O14_24216_L1"  % Verified
        xlimits = [-30 50];
        strainlimits = [-20 70];
    case "O14_24216_L2"  % Verified
        xlimits = [-30 10];
        strainlimits = [-20 80];
        locshow = 6491+([-6,-3,0,1,2,3,4,5,6])*60;
    case "O14_24216_L3"  % Verified
        xlimits = [-20 3];
        strainlimits = [-30 100];
    case "O15_1611_7882_R3"  % Verified
        xlimits = [-26 5];
        strainlimits = [-30 100];
    case "O15_1611_11589_R3"     % Verified
        strainlimits = [-30 100];
    case "O15_1611_13209_L3_1"   % Verified. Note no good strain reference is present
        xlimits = [-2 5];
    case "O15_1611_13209_L3_2"   % Verified
        xlimits = [-18 -5];
        strainlimits = [-30 80];
    case "O15_1611_13209_L4"     % Verified.
        xlimits = [-6.2 5.2];
        strainlimits = [-20 80];
        locshow = 997+([-5,-2,-1,0,1,2,3,4])*60;
    case "O15_1611_15512_L1"     % Verified
        xlimits = [-11 5];
        strainlimits = [-30 80];
        locshow = 680+(-4:1:4)*60;
    case "O15_1611_19083_L2"     % Verified
        xlimits = [-16 35];
        strainlimits = [-30 110];
    case "O15_1611_19083_L3"     % Verified
        strainlimits = [-20 110];
        xlimits = [-25 5];
        % locshow = 2935-4+(-23:3:5)*60;
        % locshow = 1279:1:3239;
    case "O15_1611_22417_L2"     % Verified
        xlimits = [-15 15];
        strainlimits = [-30 80];
    case "O15_1611_22417_L3"     % Verified
        xlimits = [-30 4];
        strainlimits = [-30 110];
        % locshow = 390:3076;
    case "O15_1611_22417_L4"
        %     at 9 seconds
        xlimits = [-10 5];
        strainlimits = [-30 80];
end

mode = "absolute";    % Options: "absolute" or "normalized"
bininterval = 0;      % Number of seconds for bininterval. Set 0 for OFF
plotpeak = false;        % Whether to plot peak on curvature plot
showpoints = false;      % Whether to show red lines for where the points are, and their labels
timeaxis = "s";     % How to set the time axis: "frame" or "s" or ???
armtouchpoint_mm = 0;   % This will be replaced if available

%% Data load

matpath = octo_InitializeData(clip, output_data_path, input_data_path);
load(matpath,'armlift','framerate','armtrimloc','curvframestart',...
    'curvframeend','armtouchpoint_mm','dropboxpath');
preprocessdata = octo_PreProcess(clip, dropboxpath, input_data_path);
load(preprocessdata,'ptorder','guidepts','numframes','ptdatmm',...
    'guideptind','arr','numpsplineseg','ptorderNG','ptnames')

% This is where the actual analysis is performed
analysisdatapath = octo_AnalyzeSegments(clip, dropboxpath, input_data_path);
load(analysisdatapath,'cumdist','curvdat','curvdatbinned','timax',...
    'segdists','segdistbig','curvPeakInd','ptdistarr','splpts');
% segdistbig = 100*segdistbig./min(segdistbig,[],1,'omitnan')-100;

% Make binned time axis
% Constrain time axis
curvframeend = min(curvframeend,numel(timax));
switch timeaxis
    case "frame"
        timax_mod = curvframestart:curvframeend;    % Use this line for frame numbers
        timax = 1:numel(timax);
        bininterval = 0;    % Binning off in this case
    case "s"
        timax_mod = timax(curvframestart:curvframeend);
end
timax_orig = timax;
segdistbig_orig = segdistbig;
if bininterval > 0
    timax_binned = (floor(min(timax)/bininterval)*bininterval):...
        bininterval:(ceil(max(timax)/bininterval)*bininterval);
    % Calculate in new grid to collect all data
    strainmat = nan(size(segdistbig,2),numel(timax_binned));
    curvmat = nan(numel(timax_binned),size(curvdatbinned,2));
    for ll = 1:size(segdistbig,2)
        % Relate input positions to sample areas
        for tt = 1:numel(timax_binned)
            tind = find(and(timax>timax_binned(tt),timax<=timax_binned(min(end,tt+1))));
            temparr = segdistbig(tind,ll);
            strainmat(ll,tt) = median(temparr(:),'omitnan');
        end
    end
    for ll = 1:size(curvdatbinned,2)
        % Relate input positions to sample areas
        for tt = 1:numel(timax_binned)
            tind = find(and(timax>timax_binned(tt),timax<=timax_binned(min(end,tt+1))));
            temparr = curvdatbinned(tind,ll);
            curvmat(tt,ll) = median(temparr(:),'omitnan');
        end
    end
    segdistbig = strainmat';
    timax = timax_binned;
    ptdistarr = ptdistarr(1:numel(timax_binned),:);
    curvdatbinned = curvmat;
else
end

% Array for drawing point locations if desired
[yyy,~] = meshgrid(timax,1:numel(ptorderNG));
disp('Processing complete. Now execute the desired section in octo_VisualizeCurvatureAndStrainIndividual.m')
return

%% Plot curvature surface

fsplcurv = figure(8);
fsplcurv.Position = [100 100 540 480];

ptnames2 = ptnames(~guideptind);

% Set limits for labeling convenience
if exist("xlimits","var")
    xl1 = xlimits(1);
    xl2 = xlimits(2);
else
    xl1 = timax_mod(1);
    xl2 = timax_mod(end);
end

switch mode
    case "normalized"
        % Normalized
        imm = imagesc(timax,cumdist(1:end-1)/max(cumdist),curvdatbinned');  % Normalized
        ylabel('Relative position along arm')
        c = colorbar;
        c.Label.String = '(Bend radius)^{-1} (mm^{-1})';
        clim([-0.2 25])
        ylim([0 1]);
        hold on
        if showpoints
            plot(yyy',ptdistarr/max(cumdist),'-r'); % Normalized
            % Add line labels
            for i=1:numel(ptorderNG)
                text(xl2-0.01*(xl2-xl1),...
                    ptdistarr(end,i)/max(cumdist),ptnames2(i),...
                    'Horizontalalignment','right','VerticalAlignment','bottom',...
                    'Color','red')  % Normalized
            end
        end
        if plotpeak
            % Plot curvature peak (normalized)
            plot(timax_orig(~isnan(curvPeakInd)),interp1(1:numel(cumdist),...
                cumdist,curvPeakInd(~isnan(curvPeakInd)))/max(cumdist),'.r');
        end

    case "absolute"
        % armtouchpoint_mm = 0;
        % Not normalized
        imm = imagesc(timax,cumdist(1:end-1)-armtouchpoint_mm,curvdatbinned');
        ylabel('Average position along arm (mm)')
        c = colorbar;
        c.Label.String = '(Bend radius)^{-1} (mm^{-1})';
        clim([-0.2 25])
        hold on
        if showpoints
            plot(yyy',ptdistarr-armtouchpoint_mm,'-r');
            % Add line labels
            for i=1:numel(ptorderNG)
                text(xl2-0.01*(xl2-xl1),...
                    ptdistarr(end,i)-armtouchpoint_mm,ptnames2(i),...
                    'Horizontalalignment','right','VerticalAlignment','bottom',...
                    'Color','red')  % Not normalized
            end
        end
        if plotpeak
            % Plot curvature peak (not normalized)
            plot(timax_orig(~isnan(curvPeakInd)),interp1(1:numel(cumdist),...
                cumdist,curvPeakInd(~isnan(curvPeakInd)))-armtouchpoint_mm,'.r');
        end
end
switch timeaxis
    case "frame"
        xlabel('Frame number')
    case "s"
        xlabel('Time from substrate release (s)');
end

hold off
xlim([timax_mod(1) timax_mod(end)])
if exist("xlimits","var")
    xlim(xlimits);
end
% title(clip,'Interpreter','none')
colormap parula
cc=colormap;
cc(1,:) = [0.9 0.9 0.9];
colormap(cc)
curvax = gca;
curvax.FontSize = 20;
curvax.YDir = "normal";
curvax.Position(3)=curvax.Position(3)*0.95;

%% Plot shift in curvature peak over time

curvshiftfig = figure(11);
curvshiftfig.Position = [550 100 440 480];
plot(timax_orig(~isnan(curvPeakInd)),interp1(1:numel(cumdist),...
    cumdist,curvPeakInd(~isnan(curvPeakInd)))-armtouchpoint_mm,'.r');
ylim([-80 20])
xlabel('Time from substrate release (s)')
ylabel('Average position along arm (mm)')
title(clip,'Interpreter','none')
optimizeFig;

%% Show stretching with similar graph as curvature surface above

fstretch2D = figure(9);
fstretch2D.Position = [100 680 440 480];

% NOTE: SELECT YOUR FAVORITE REFERENCE HERE!
% (1) Use minimum as reference
% segdistbig = 100*segdistbig./min(segdistbig,[],1,'omitnan')-100;
% cbarstring = "Strain (% from minimum)";

% (2) Use time 0 for strain reference
% tind0 = armlift;
% % If not all segments are populated, throw warning.
% temp = ~isnan(sum(segdists,2));
% if temp(tind0) == 0
%     warning('Not all segments populated');
% end
% segdistbig = 100*segdistbig./segdistbig(tind0,:);
% cbarstring = "Strain (% from t=0)";

% (3) Use average in window -20 to -9 for strain reference
tref = (timax>-20) & (timax<=-9);
segdistbig_orig = 100*segdistbig_orig./median(segdistbig(tref,:),1,'omitnan')-100;
segdistbig = 100*segdistbig./median(segdistbig(tref,:),1,'omitnan')-100;
cbarstring = "Strain (% from median at -20<t<=-9)";

% % (3b) Use average in window -10 to -5 for strain reference
% tref = (timax>-10) & (timax<=-5);
% segdistbig_orig = 100*segdistbig_orig./median(segdistbig(tref,:),1,'omitnan')-100;
% segdistbig = 100*segdistbig./median(segdistbig(tref,:),1,'omitnan')-100;
% cbarstring = "Strain (% from median at -10<t<=-5)";

% % (4) Use average in window 0 to 2 for strain reference
% tref = (timax>0) & (timax<=2);
% segdistbig_orig = 100*segdistbig_orig./median(segdistbig(tref,:),1,'omitnan')-100;
% segdistbig = 100*segdistbig./median(segdistbig(tref,:),1,'omitnan')-100;
% cbarstring = "Strain (% from median at 0<t<2)";

switch mode
    case "normalized"
        imm = imagesc(timax,ptdistarr(1,:)/max(cumdist),segdistbig');     % Normalized

        ylabel('Relative position along arm')
        c = colorbar;

        hold on
        if showpoints
            plot(yyy',ptdistarr/max(cumdist),'-r'); % Normalized
            % Add line labels
            for i=1:numel(ptorderNG)
                text(xl2-0.01*(xl2-xl1),...
                    ptdistarr(end,i)/max(cumdist),ptnames2(i),...
                    'Horizontalalignment','right','VerticalAlignment','bottom',...
                    'Color','red')  % Normalized
            end
        end
    case "absolute"
        imm = imagesc(timax,cumdist(1:end-1)-armtouchpoint_mm,segdistbig');     % Not normalized
        ylabel('Average position along arm (mm)')
        c = colorbar;

        hold on
        if showpoints
            plot(yyy',ptdistarr-armtouchpoint_mm,'-r');
            % Add line labels
            for i=1:numel(ptorderNG)
                text(xl2-0.01*(xl2-xl1),...
                    ptdistarr(end,i)-armtouchpoint_mm,ptnames2(i),...
                    'Horizontalalignment','right','VerticalAlignment','bottom',...
                    'Color','red')  % Not normalized
            end
        end
end
switch timeaxis
    case "frame"
        xlabel('Frame number')
    case "s"
        xlabel('Time from substrate release (s)');
end

colormap(cc);
c.Label.String = cbarstring;
if exist("strainlimits","var")
    clim(strainlimits);
else
    strainlimits = clim;
end
hold off
xlim([timax_mod(1) timax_mod(end)])
if exist("xlimits","var")
    xlim(xlimits);
end

% Fix limits so that only NaN is displayed dark
cl1 = c.Limits;
c.Limits = cl1;
imm.CData(imm.CData<=(cl1(1))) = cl1(1)+diff(cl1)/size(cc,1);
imm.CData(isnan(imm.CData)) = cl1(1);

% title(clip,'Interpreter','none')
stretchax = gca;
stretchax.FontSize = 20;
stretchax.YDir = "normal";
stretchax.Position(3)=stretchax.Position(3)*0.95;

% % Output stretch and curvature data
% outfile = ['temp_MATLAB' filesep char(clip) '_stretchcurve.mat'];
% save(outfile,'stretchax','curvax','cumdist',...
%     'timax','curvdat','ptdistarr','segdistbig');


%% Box plot
segstretchfig = figure(4);
segdists2 = movmean(segdists,50,1);
boxplot(100*segdists2./min(segdists2))
xlabel('Segment number')
ylabel('Segment stretch statistics (%)')
% title('Box plot of segment stretch (relative to first frame)')
xticks(1:numel(ptorder));
xticklabels(cellstr(num2str(ptorder')));
segstretchfig.Position = [550 680 640 400];
ylim([0 600])
optimizeFig;


%% Show arm curves, colored by TIME

% numcurves = 20;     % Number of curves to show
% locshow = unique(round(linspace(1,curvframeend-curvframestart+1,numcurves)));

if ~exist("locshow","var")
    locshow = 1:180:curvframeend-curvframestart+1;      % This worked if not binned
end

% numclrs = numel(locshow);
numclrs = (max(locshow)-min(locshow))/60+1;
clrs = turbo(numclrs+2);
clrs = clrs(2:end-1,:); % First and last color are bad

farmc = figure(10);
farmc.Position = [600 720 640 500];
splptstrim = splpts(:,:,curvframestart:curvframeend);
ptdatmmtrim = ptdatmm(:,:,curvframestart:curvframeend);
ref = permute(splptstrim(1,:,locshow),[3,2,1]);
% ref = ref*0;  % This line forces "world frame"
numsplpts = size(splpts,1);
locshowvis = [];
for i=1:numel(locshow)
    locshowvis(i) = false;
    for j=1:numsplpts
        % Draw point by point
        if ~isnan(splptstrim(j,1,locshow(i)))
%             thisclr = clrs(i,:);
            thisclr = clrs(round((locshow(i)-min(locshow))/60+1),:);
            plot3(splptstrim(j,1,locshow(i))-ref(i,1),...
                -(splptstrim(j,3,locshow(i))-ref(i,3)),...
                splptstrim(j,2,locshow(i))-ref(i,2),...
                '.','MarkerSize',8,'Color',thisclr)

            hold on
            locshowvis(i) = true;
        else
            
        end
    end
    % % plot3(splpts(:,1,locshow(i))-ref(i,1),splpts(:,3,locshow(i))-ref(i,3),splpts(:,2,locshow(i))-ref(i,2),'-k','linewidth',2,'Color',clrs(i,:))
    % % hold on
    plot3(ptdatmmtrim(:,1,locshow(i))-ref(i,1),...
        -(ptdatmmtrim(:,3,locshow(i))-ref(i,3)),...
        ptdatmmtrim(:,2,locshow(i))-ref(i,2),'.r','MarkerSize',14)
%     title(num2str((locshow(i)-armlift)/60))
%     pause(1)
end
hold off

xlabel('x (mm)')
ylabel('z (mm)')
zlabel('y (mm)')

axis equal
% title(['Interval between curves: ' num2str(median(diff(locshow))/60) ' seconds'])

% surf(300+(1:size(segdistbig,1)),300+(1:size(segdistbig,2)),segdistbig);
colormap(clrs)
c = colorbar;
clim([(locshow(1)-armlift)/60-0.5 (locshow(end)-armlift)/60+0.5])
c.Label.String = 'Time from release (s)';
c.Ticks = round((locshow(locshowvis>0)-armlift)/60);

% title(['Sample arm poses from clip ' clip],'Interpreter','none');
ax2 = gca;
ax2.XGrid = "on";
ax2.YGrid = "on";
ax2.ZGrid = "on";
ax2.FontSize = 20;
view([-30,10]);
farmc.Position(3) = 756;
farmc.Position(4) = 509;
xl = xlim; yl = ylim; zl = zlim;

%% Show arm curves, colored by stretch

% numcurves = 20;     % Number of curves to show
% locshow = unique(round(linspace(1,curvframeend-curvframestart+1,numcurves)));

if ~exist("locshow","var")
    locshow = 1:180:curvframeend-curvframestart+1;      % This worked if not binned
    locshow = 1:30:curvframeend-curvframestart+1;      % This worked if not binned
    locshowold = locshow;
else
    locshowold = locshow;
    locshow = curvframestart:30:curvframeend-curvframestart+1;
end
% locshow = 920:10:1150;
% locshow = 2600:10:2800;
locshow = 390:180:3076;

numclrs = 100;
clrs = parula(numclrs);

farmc = figure(10);
pos11 = [600 720 640 600];
farmc.Position = pos11;
splptstrim = splpts(:,:,curvframestart:curvframeend);
ptdatmmtrim = ptdatmm(:,:,curvframestart:curvframeend);
ref = permute(splptstrim(1,:,locshow),[3,2,1]);
% ref = ref*0;  % This line forces "world frame"
numsplpts = size(splpts,1);
for i=1:numel(locshow)
    for j=1:numsplpts
        % Draw point by point
        if ~isnan(splptstrim(j,1,locshow(i)))
            strainval = segdistbig_orig(curvframestart+locshow(i)-1,j);
            if strainval>strainlimits(2)
                strainval = strainlimits(2);
            elseif strainval<strainlimits(1)
                strainval = strainlimits(1);
            end
            if isnan(strainval)
                thisclr = [0.5 0.5 0.5];
            else
                thisclr = interp1(linspace(strainlimits(1),strainlimits(2),numclrs),clrs,strainval);
            end
            plot3(splptstrim(j,1,locshow(i))-ref(i,1),...
                -(splptstrim(j,3,locshow(i))-ref(i,3)),...
                splptstrim(j,2,locshow(i))-ref(i,2),...
                '.','MarkerSize',5,'Color',thisclr)

            hold on
        end
    end
    % % plot3(splpts(:,1,locshow(i))-ref(i,1),splpts(:,3,locshow(i))-ref(i,3),splpts(:,2,locshow(i))-ref(i,2),'-k','linewidth',2,'Color',clrs(i,:))
    % % hold on
    plot3(ptdatmmtrim(:,1,locshow(i))-ref(i,1),...
        -(ptdatmmtrim(:,3,locshow(i))-ref(i,3)),...
        ptdatmmtrim(:,2,locshow(i))-ref(i,2),'.r','MarkerSize',12)
%     title(num2str((locshow(i)-armlift)/60))
%     pause(1)
end
hold off

xlabel('x (mm)')
ylabel('z (mm)')
zlabel('y (mm)')

axis equal
% title(['Interval between curves: ' num2str(median(diff(locshow))/60) ' seconds'])

% surf(300+(1:size(segdistbig,1)),300+(1:size(segdistbig,2)),segdistbig);
colormap(cc)
c = colorbar;
clim(strainlimits)
c.Label.String = cbarstring;

% title(['Sample arm poses from clip ' clip],'Interpreter','none');
ax2 = gca;
ax2.XGrid = "on";
ax2.YGrid = "on";
ax2.ZGrid = "on";
ax2.FontSize = 20;
if clip == "O15_1611_22417_L3"
    view([-35,8]);   % Use for 22417_L3
    ax2.YLabel.Position(1) = -40;
    ax2.YLabel.Position(2) = -60;
    xlim([-5 105])
%     ylim([-36,24]);
    ylim([-45 26]);
    zlim([-61 5]);
else
    view([35,8]);   % Use for 19083_L3
    zl = zlim+[-10 10]; zlim(zl);
end
% farmc.Position(3) = pos11(3);
% farmc.Position(4) = pos11(4);
xl = xlim; yl = ylim; zl = zlim;

% locshow = locshowold;

%% Show arm curves, colored by stretch - supp mat video - be sure to run section above first

if ~exist("locshow","var")
    locshow = 1:180:curvframeend-curvframestart+1;      % This worked if not binned
end

% locshow = 1279:1:3239;    % 19083
locshow = 390:3076;     % 22417_L3
% locshow = curvframestart:180:curvframeend-curvframestart+1;

numclrs = 100;
clrs = parula(numclrs);

farmc = figure(10);
farmc.Position = pos11;
splptstrim = splpts(:,:,curvframestart:curvframeend);
ptdatmmtrim = ptdatmm(:,:,curvframestart:curvframeend);
ref = permute(splptstrim(1,:,locshow),[3,2,1]);
ref = fixref(ref);

% ref = ref*0;  % This line forces "world frame"
numsplpts = size(splpts,1);

outvid3d = [dropboxpath 'visualizations' filesep char(clip) '_suppvid_strain.mp4'];
vw = VideoWriter(outvid3d,'MPEG-4');
open(vw);

% Initialize view
% farmc = figure(10);
axis(ax2,'manual')
thisfig=gcf;
hold(ax2,"off")
k=1;
for i=1:numel(locshow)
    if i==1
        for j=1:numsplpts
            p(j)=plot3(ax2,nan,nan,nan,'.','MarkerSize',5,'Color',[0,0,0]);
            hold on
        end
    end
    for j=1:numsplpts
        % Draw point by point
        if ~isnan(splptstrim(j,1,locshow(i)))
            strainval = segdistbig_orig(curvframestart+locshow(i)-1,j);
            if strainval>strainlimits(2)
                strainval = strainlimits(2);
            elseif strainval<strainlimits(1)
                strainval = strainlimits(1);
            end
            if isnan(strainval)
                thisclr = [0.5 0.5 0.5];
            else
                thisclr = interp1(linspace(strainlimits(1),strainlimits(2),numclrs),clrs,strainval);
            end
            if i==1
            p(j)=plot3(ax2,splptstrim(j,1,locshow(i))-ref(i,1),...
                -(splptstrim(j,3,locshow(i))-ref(i,3)),...
                splptstrim(j,2,locshow(i))-ref(i,2),...
                '.','MarkerSize',5,'Color',thisclr);
%             plot3(ax2,splptstrim(j,1,locshow(i))-ref(i,1),...
%                 -(splptstrim(j,3,locshow(i))-ref(i,3)),...
%                 splptstrim(j,2,locshow(i))-ref(i,2),...
%                 '.','MarkerSize',5,'Color',thisclr);
            hold on
            else
                p(j).XData = splptstrim(j,1,locshow(i))-ref(i,1);
                p(j).YData = -(splptstrim(j,3,locshow(i))-ref(i,3));
                p(j).ZData = splptstrim(j,2,locshow(i))-ref(i,2);
                p(j).Color = thisclr;
            end
        else
            p(j).XData = nan;
            p(j).YData = nan;
            p(j).ZData = nan;
        end
    end
    % plot3(splpts(:,1,locshow(i))-ref(i,1),splpts(:,3,locshow(i))-ref(i,3),splpts(:,2,locshow(i))-ref(i,2),'-k','linewidth',2,'Color',clrs(i,:))
    
    if i==1
        p(j+1) = plot3(ax2,ptdatmmtrim(:,1,locshow(i))-ref(i,1),...
            -(ptdatmmtrim(:,3,locshow(i))-ref(i,3)),...
            ptdatmmtrim(:,2,locshow(i))-ref(i,2),'.r','MarkerSize',12);
        hold off
    else
        p(j+1).XData = ptdatmmtrim(:,1,locshow(i))-ref(i,1);
        p(j+1).YData = -(ptdatmmtrim(:,3,locshow(i))-ref(i,3));
        p(j+1).ZData = ptdatmmtrim(:,2,locshow(i))-ref(i,2);
    end

%     plot3(ax2,ptdatmmtrim(:,1,locshow(i))-ref(i,1),...
%         -(ptdatmmtrim(:,3,locshow(i))-ref(i,3)),...
%         ptdatmmtrim(:,2,locshow(i))-ref(i,2),'.r','MarkerSize',12);
%     hold off

    % Update view settings on first round
    if i==1
        xlabel('x (mm)')
        ylabel('z (mm)')
        zlabel('y (mm)')

        axis equal
        % title(['Interval between curves: ' num2str(median(diff(locshow))/60) ' seconds'])

        % surf(300+(1:size(segdistbig,1)),300+(1:size(segdistbig,2)),segdistbig);
        colormap parula
        c = colorbar;
        clim(strainlimits)
        c.Label.String = cbarstring;

        % title(['Sample arm poses from clip ' clip],'Interpreter','none');
        ax2 = gca;
        ax2.XGrid = "on";
        ax2.YGrid = "on";
        ax2.ZGrid = "on";
        ax2.FontSize = 20;
        temp = xlim;
        xl(1) = min(xl(1),temp(1));
        xl(2) = max(xl(2),temp(2));
        xl
        temp = ylim;
        yl(1) = min(yl(1),temp(1));
        yl(2) = max(yl(2),temp(2));
        temp = zlim;
        zl(1) = min(zl(1),temp(1));
        zl(2) = max(zl(2),temp(2));
        axis manual
        xlim(ax2,xl+[-5,5]);
        ylim(ax2,yl+[-5,5]);
        zlim(ax2,zl+[-10,5]);
        if clip == "O15_1611_22417_L3"
            view([-35,8]);   % Use for 22417_L3
            ax2.YLabel.Position(1) = -20;
            ax2.YLabel.Position(2) = -40;
        else
            view([35,8]);   % Use for 19083_L3
        end
    end

    drawnow;
%     pause(1);
    % Grab overlaid frame
    frame = getframe(thisfig);
    frame = imresize(frame.cdata,2);
%     frame=frame.cdata;
    i
%     size(frame)
    writeVideo(vw,frame);
end
close(vw)

%% Functions
function outref = fixref(ref)
% Interpolate values in reference point as needed
if isnan(ref(1,1))
    ind = find(~isnan(ref(:,1)),1,'first');
    ref(1:ind,:) = repmat(ref(ind,:),ind,1);
end
if isnan(ref(end,1))
    ind = find(~isnan(ref(:,1)),1,'last');
    ref(ind:end,:) = repmat(ref(ind,:),size(ref,1)-ind+1,1);
end
ind = isnan(ref(:,1));
indx = 1:size(ref,1);
ref(ind,1) = interp1(indx(~ind),ref(~ind,1),find(ind));
ref(ind,2) = interp1(indx(~ind),ref(~ind,2),find(ind));
ref(ind,3) = interp1(indx(~ind),ref(~ind,3),find(ind));
outref = ref;
end
