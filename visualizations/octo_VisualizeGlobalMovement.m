%% Plot temporal "global" movement data derived from multiple camera angles
% Uses data from Google Docs to make plots that can help detect patterns

%% Inputs
clear all
gsheetID = ['1BfUr0cvolhZkh_jn-tkxuUsW' ...
    'mn_M_41eYFIOi4uo-9M'];
basetime = datetime('20220826T231332','InputFormat','yyyyMMdd''T''HHmmss');
framerate = 29.97;

%% Load Google Sheet data

localFile = 'gsheetdatacache.mat';
try
    sheetdata = getGoogleSheetData(gsheetID);
    % Save data in a temporary file for offline access
    save(localFile,'sheetdata');
catch
    % Load local cached file if Google Sheet is not available
    load(localFile,'sheetdata');
    disp('Local data loaded');
end
tabl = sheetdata(:,1:end-1);

% Re-order if desired
% reord = [1 2 9 10 3 4 11 12 5 6 13 14 7 8 15 16];   % L1 R1 L2 R2
reord = [1:8 15 16 13 14 11 12 9 10];   % L1 L2 ... L4 R4 R3
tabl = tabl(:,reord);

for i=1:2:size(tabl,2)
  tabl.(i) = str2double(tabl{:,i});
end
tablarr=table2array(tabl(:,1:2:end));


%% Generate full temporal data

startfr = min(tablarr(:));
tarray = startfr:max(tablarr(:));
numarms = size(tablarr,2);
globarm = nan(numel(tarray),numarms);
% 1 is down, 0 is up
for i = 1:size(tablarr,2)
    for j=1:find(~isnan(tablarr(:,i)),1,'last')
        switch tabl{j,i*2}{:}
            case 'D'
                val = 1;
            case 'U'
                val = 0;
            otherwise
                val = nan;
        end
        startind = tablarr(j,i)-startfr+1;
        if j==find(~isnan(tablarr(:,i)),1,'last')
            endind = numel(tarray);
        else
            endind = tablarr(j+1,i)-startfr+1;
        end
        globarm(startind:endind,i) = val;
    end
end

realtimes = basetime + seconds(tarray/framerate);


%% Plot

f96 = figure(96);
f96.Position = [160 800 1160 400];
timax = seconds((realtimes)-(basetime));
timeref = timax(1);
imagesc(timax-timeref,1:numarms,globarm')
% imagesc(1:numel(realtimes),1:numarms,globarm')
caxis([-1 1]);
touchcolor = [0 0 0];
notouchcolor = [1 1 1];
nodatacolor = [0.85 0.85 0.85];
cc = [nodatacolor;notouchcolor;touchcolor];
colormap(cc)
% xlabel(['Time (seconds from ' datestr(basetime) ')'])
xlabel('Time (s)')


% armnames = sheetdata.Properties.VariableNames(1:2:end-1);
armnames = tabl.Properties.VariableNames(1:2:end);
for i=1:numel(armnames)
    armnames(i) = {armnames{i}(1:2)};
end
yticklabels(armnames)
% title(['Black = arm touching substrate, White = not touching substrate, ' ...
%     'Grey = no data, Red = EyeRIS data analyzed']);

% Add markings for available data
clips = {"octo_1611_7318_L1","octo_1611_7318_L2",...
    "octo_1611_7882_R1","octo_1611_7882_R3","octo_1611_7882_R4",...
    "octo_1611_11589_R3","octo_1611_11589_R4",...
    "octo_1611_13209_L2","octo_1611_13209_L3",...
    "octo_1611_15512_L1",...
    "octo_1611_19083_L1","octo_1611_19083_L2","octo_1611_19083_L3",...
    "octo_1611_22417_L1","octo_1611_22417_L2","octo_1611_22417_L3","octo_1611_22417_L4"};
for i=1:numel(clips)
    % Locate valid EyeRIS data for each clip
    clip = clips{i};
    corrdata = octo_DataCorrection(clip);
    load(corrdata,'ptdatmm','quadstartframe');
    hold on
    temp = squeeze(sum(~isnan(ptdatmm(:,3,:)),1))>1;    % Indices of times with at least two points with z data
    clipstart = quadstartframe/29.97;
    tempnm = char(clip);
    temp2=find(temp);       % Locations of desired times
    temp3=[1 find(diff(find(temp))>1)'+1 numel(temp2)];
    for j = 1:(numel(temp3)-1)
        startfr = clipstart+(temp2(temp3(j))-1)/60-timeref;
        endfr = clipstart+(temp2(temp3(j+1)-1)-1)/60-timeref;
        ypos = find(strcmp(armnames,tempnm(end-1:end)));
        rectangle("Position",[startfr, ypos-0.5, endfr-startfr,1],...
            'EdgeColor','r',"LineWidth",2);
    end
end
hold off
xl = xlim;
xlim([0 xl(2)])
% optimizeFig
thisax = gca;
thisax.FontSize = 19;
thisax.Clipping = 'off';
hold on
walkstartfr = 11750;
walkstart = walkstartfr/framerate-timeref;
plot([walkstart walkstart],[0.5 9.7],'--','Color',[0.5 0.5 0.5])
text(walkstart,9.7,'--> crawling','Color',[0.5 0.5 0.5],...
    'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',17)
% text(walkstart,9.5,'no locomotion <--','Color',[0.5 0.5 0.5],...
%     'HorizontalAlignment','right','VerticalAlignment','baseline','FontSize',12)
hold off
dfdf
%% Simple statistics

intv = timax(2)-timax(1);
for i = 1:8
    currarm = globarm(walkstartfr-tarray(1)+1:end,i);   % Use only section where walking

    % Reduce array to start of the first down to start of last down
    temp = diff(currarm);
    ind1 = find(temp==1,1,'first')+1;
    ind2 = find(temp==1,1,'last');
    if ind1>ind2
        warning(['Problem in ' armnames{i}]);
        currarm = nan;
    else
    currarm = currarm(ind1:ind2);
    end
    
    % Get percentage of time the arm is down
    downperc(i) = 100*sum(currarm==1)/sum(~isnan(currarm));
    downpercmax(i) = 100*(sum(currarm==1)+sum(isnan(currarm)))/numel(currarm);
    downpercmin(i) = 100*(sum(currarm==1))/numel(currarm);

%     % Make temporary array that fills NaNs with most recent finite value
%     temp = currarm;
%     for j=2:find(~isnan(temp),1,'last')
%         if isnan(temp(j))
%             temp(j) = temp(j-1);
%         end
%     end

    % Make temporary array that fills NaNs, through interpolation
    % (Switching halfway if transition occurs)
    temp = currarm;
    nanind = isnan(temp);
    if sum(nanind)>0 && sum(~nanind)>2
        tempind = 1:numel(temp);
        temp(nanind) = round(interp1(tempind(~nanind),temp(~nanind),tempind(nanind)));
    end

    % Get frequency of substrate release (# of releases / total tracked
    % time)
    numups(i) = sum(diff(temp)==-1);
    release_frequency(i) = numups(i) / ...
        (numel(currarm)*intv);

    downtimes = [];
    cycleduration = [];
    % Label downtimes
    temp2 = diff(temp); temp2(temp2<0) = 0;
    temp2 = cumsum([1; temp2]);
    temp2(temp<1) = NaN;
    for j = 1:numups(i)-1
        downtimes(j) = sum(temp2==j)*intv;
        cycleduration(j) = (find(temp2==j+1,1,'first')-find(temp2==j,1,'first'))*intv;
    end
    j=j+1;
    downtimes(j) = sum(temp2==j)*intv;
    cycleduration(j) = (numel(temp2)+1-find(temp2==j,1,'first'))*intv;
    cycledurationmean(i) = mean(cycleduration);
    cycledurationstd(i) = std(cycleduration);
    avdowntime(i) = mean(downtimes); 
    avdowntimestd(i) = std(downtimes);

    downperc(i) = 100*median(downtimes./cycleduration);


    % Calculate duration of downtime (stats), both optimistic/max (nans =
    % down) and conservative/min (nans = up). % Do not count down times at
    % beginning and end of track

end
outtable = table(downperc',downpercmin',downpercmax',numups',release_frequency',avdowntime',avdowntimestd');
outarr = table2array(outtable);


%% Add some statistics to plot

% Move graph to left to make room for stats
thisax.Position(1) = 0.05;
thisax.Position(3) = 0.70;

% Stats titles
title1x = max(timax)-timeref+20;
t1x = text(title1x,0,'Cycle duration (s)','FontSize',14);
title2x = t1x.Position(1)+ 140;
t2x = text(title2x,0,'Down time (%)','FontSize',14);

% Table lines
line([title2x-20 title2x-20],[0 8.5],'Color',[0.5 0.5 0.5])
for i=0:8
    line([max(timax)-timeref title2x+90],0.5+[i i],'Color',[0.5 0.5 0.5])
end

% Add stat data
for i=1:8
    text(title1x+20,i,[num2str(cycledurationmean(i),'%.1f') ' (' num2str(cycledurationstd(i),'%.1f') ')'])
    text(title2x+30,i,num2str(downperc(i),'%.1f'))
end

% title('Need to verify how unknowns are treated in stats!')

sdfsdfs
%% Spectral analysis

f95 = figure(95);
f95.Name = 'Fourier transforms';

Fs = 1/(timax(2)-timax(1)); % Sample frequency
L = 1e5; % Signal length (will pad zeros if needed)
f = Fs*(0:(L/2))/L; % Frequency array
P1all = zeros(L/2+1,numarms);

for i=1:8

% Do Fourier transform
currarm = globarm(:,i);
currarm(isnan(currarm)) = 0;
Y = fft(currarm,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1all(:,i) = P1;

% Plot
ordr = [1:4 8:-1:5];
subplot(2,4,ordr(i));
plot(f,P1);
title(['Arm ' armnames{i}])
xlabel("f (Hz)")
ylabel("|P1(f)|")
xlim([0 0.1])
end


% f94 = figure(94);
% f94.Name = 'Walsh-Hadamard transform';
% 
% Fs = 1/(timax(2)-timax(1)); % Sample frequency
% L = 1e5; % Signal length (will pad zeros if needed)
% f = Fs*(0:(L/2))/L; % Frequency array
% P1all = zeros(L/2+1,numarms);
% 
% for i=1:8
% 
% % Do Fourier transform
% currarm = globarm(:,i);
% currarm(isnan(currarm)) = 0;
% y = fwht(currarm);
% 
% % Plot
% ordr = [1:4 8:-1:5];
% subplot(2,4,ordr(i));
% plot(abs(y))
% xlabel('Sequency index')
% ylabel('Magnitude')
% title(['Arm ' armnames{i}])
% % xlabel("f (Hz)")
% % ylabel("|P1(f)|")
% % xlim([0 0.1])
% end

