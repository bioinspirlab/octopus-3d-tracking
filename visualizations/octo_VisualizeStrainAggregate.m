%% Information

% Create visualization that looks at the peaks in curvature along each arm,
% and plots them all together

clear all
close all
%% Inputs

% clips = {"O10_6946_L1",...
%     "O10_7875_R1",...
%     "O14_24216_L1",...
%     "O14_24216_L2",...
%     "O14_24216_L3",...
%     "O15_1611_7882_R3",...
%     "O15_1611_11589_R3",...
%     "O15_1611_13209_L3_1",...
%     "O15_1611_13209_L3_2",...
%     "O15_1611_13209_L4",...
%     "O15_1611_15512_L1",...
%     "O15_1611_19083_L2",...
%     "O15_1611_19083_L3",...
%     "O15_1611_22417_L2",...
%     "O15_1611_22417_L3",...
%     "O15_1611_22417_L4"};

% % Arm 3 and 4 only
% clips = {"O14_24216_L3",...
%     "O15_1611_7882_R3",...
%     "O15_1611_11589_R3",...
%     "O15_1611_13209_L3_1",...
%     "O15_1611_13209_L3_2",...
%     "O15_1611_13209_L4",...
%     "O15_1611_19083_L3",...
%     "O15_1611_22417_L3",...
%     "O15_1611_22417_L4"};
% strainlimits = [-15 50];

% % Arm 2 only
% clips = {"O14_24216_L2",...
%     "O15_1611_19083_L2",...
%     "O15_1611_22417_L2"};
% % strainlimits = [-15 70];

% Arm 1 and 2
clips = {"O10_6946_L1",...
    "O10_7875_R1",...
    "O14_24216_L1",...
    "O14_24216_L2",...
    "O15_1611_15512_L1",...
    "O15_1611_19083_L2",...
    "O15_1611_22417_L2"};
strainlimits = [-15 50];

% strain_thresh = 50;     % Percentage strain above which "high strain" is defined
mode = "scaled";    % Y axis options: "absolute" (mm) or "scaled" (by head width)
x_mode = "scaled";  % X axis options: "absolute" (s) or "scaled" (by gait duration)

% Set axis intervals
switch mode
    case "absolute"
        posarr = -60:5:20;      % Use for "absolute"
    case "scaled"
        posarr = -1.1:0.1:0.4;  % Use for "scaled"
end
switch x_mode
    case "absolute"
        bininterval = 2;        % Number of seconds for bininterval. Set 0 for OFF
        timarr = -28:bininterval:4;
    case "scaled"
        bininterval = 0.05;        % Size of bininterval. Set 0 for OFF
        timarr = -0.7:bininterval:0.1;
        timarr = -0.5:bininterval:0.4;  % Used for arms 1 and 2!
end


%% Processing

clrs = jet(numel(clips));
strainmat = nan(numel(posarr),numel(timarr),numel(clips));

for i = 1:numel(clips)
%     for i=2
    matpath = octo_InitializeData(clips{i});
    load(matpath,'armlift','armtrimloc','curvframestart',...
        'curvframeend','armtouchpoint_mm','headwidth','touchduration',...
        'gaitduration');
    
    % This is where the actual analysis is performed
    analysisdatapath = octo_AnalyzeSegments(clips{i});
    load(analysisdatapath,'cumdist','curvdat','curvdatbinned','timax',...
        'segdists','segdistbig','curvPeakInd','ptdistarr');

    %% Plot curvature surface

    fsplcurv = figure(10);
    fsplcurv.Position = [100 100 560 440];

    if numel(curvPeakInd)~=numel(timax)
        warning('Different array lengths');
    end

    % Make x axis
%     temp = find(~isnan(curvPeakInd));
    temp = 1:numel(timax);
    temp(temp<curvframestart) = [];
    temp(temp>curvframeend) = [];

    timax_mod = timax(temp);

% % Find time 0 for curvature reference
% tind0 = find(timax>=0,1,'first');
% % If not all segments are populated, find the closest one that is.
% temp = ~isnan(sum(segdists,2));
% if temp(tind0) == 0
%     tind0 = find(temp(1:tind0)==1,1,'last');
% end
% segdistbig = segdistbig./segdistbig(tind0,:);
% % segdistbig = segdistbig./min(segdistbig,[],1,'omitnan');

% (3) Use average in window -20 to -10 for strain reference
switch x_mode
    case "absolute"
        tref = (timax>-20) & (timax<=-10);
        cbarstring = "Strain (% from median at -20<t<=-10)";
    case "scaled"
        timax = timax/gaitduration;
        tref = (timax>-0.5) & (timax<=-0.25);
        cbarstring = "Strain (% from median at -0.5<t<=-0.25)";
end
segdistbig = 100*segdistbig./median(segdistbig(tref,:),1,'omitnan')-100;

% % (3b) Use average in window -20 to -8 for strain reference
% tref = (timax>-20) & (timax<=-8);
% segdistbig = 100*segdistbig./median(segdistbig(tref,:),1,'omitnan')-100;
% cbarstring = "Strain (% from median at -20<t<-8)";

% Calculate in new grid to collect all data
    switch mode
        case "scaled"
            inppos = ((1:size(segdistbig,2))*ptdistarr(1,end)/size(segdistbig,2)-armtouchpoint_mm)/headwidth;
        case "absolute"
            inppos = (1:size(segdistbig,2))*ptdistarr(1,end)/size(segdistbig,2)-armtouchpoint_mm;
    end
    for ll = 1:numel(posarr)
        % Relate input positions to sample areas
        llind = find(and(inppos>posarr(ll),inppos<=posarr(min(end,ll+1))));
        for tt = 1:numel(timarr)
            tind = find(and(timax>timarr(tt),timax<=timarr(min(end,tt+1))));
            temparr = segdistbig(tind,llind);
            strainmat(ll,tt,i) = median(temparr(:),'omitnan');
        end
    end

imagesc(timarr,posarr,strainmat(:,:,i))
title(clips(i),'Interpreter','none')
%     hold on
pause(0)
end

%% Make aggregate plot and adjust appearance

imm = imagesc(timarr+bininterval/2,posarr+median(diff(posarr))/2,median(strainmat,3,'omitnan'));
hold off

if exist("strainlimits","var")
    clim(strainlimits);
end
ylim([min(posarr) max(posarr)])
xlim([min(timarr) max(timarr)])
ax = gca;
ax.YDir='normal';
switch x_mode
    case "scaled"
        xlabel('Time from substrate release (normalized)')
    case "absolute"
        xlabel('Time from substrate release (s)')
end
switch mode
    case "scaled"
        ylabel('Distance from attachment point (normalized)')
    case "absolute"
        ylabel('Average position along arm (mm)')
end
if exist("titlestring","var")
    title(titlestring,'Interpreter','none')
end
% legend(clips,'Interpreter','none')
% optimizeFig;
c = colorbar;
c.Label.String = cbarstring;
ax.FontSize = 15;

cmap = ax.Colormap;
cmap(1,:)=[0.9 0.9 0.9];
ax.Colormap = cmap;

% Fix limits so that only NaN is displayed dark
cl1 = c.Limits;
c.Limits = cl1;
imm.CData(imm.CData<=(cl1(1))) = cl1(1)+diff(cl1)/size(cmap,1);
imm.CData(isnan(imm.CData)) = cl1(1);

%% Also indicate how many sequences went into each data point


fdata = figure(11);
fdata.Position = [100 700 560 440];

numreplicates = sum(~isnan(strainmat),3);
imm = imagesc(timarr+bininterval/2,posarr+median(diff(posarr))/2,numreplicates);
hold off

ylim([min(posarr) max(posarr)])
xlim([min(timarr) max(timarr)])
ax = gca;
ax.YDir='normal';
switch x_mode
    case "scaled"
        xlabel('Time from substrate release (normalized)')
    case "absolute"
        xlabel('Time from substrate release (s)')
end
switch mode
    case "scaled"
        ylabel('Distance from attachment point (normalized)')
    case "absolute"
        ylabel('Average position along arm (mm)')
end
% title('Number of sequences with data','Interpreter','none')
imm.CData(imm.CData==max(numreplicates(:))) = max(numreplicates(:))+0.5;
imm.CData(imm.CData==0) = -0.5;
ax.FontSize = 15;
c = colorbar;
c.Label.String = 'Number of sequences';
c.Ticks = 0:1:max(numreplicates(:));
% c.Limits = [-0.5 max(numreplicates(:))+0.5];

% cmap = ax.Colormap;
cmap = parula(max(numreplicates(:))+1);
cmap(1,:)=[0.9 0.9 0.9];
ax.Colormap = cmap;
