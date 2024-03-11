%% Information

% Analysis per segment (after running preprocessing steps)
% MORE DETAILS

function analysisdatapath = octo_AnalyzeSegments(clip)

% Load initialization & preprocess data
matpath = octo_InitializeData(clip);
load(matpath,'dropboxpath','armlift','framerate','armtrimloc');
preprocessdata = octo_PreProcess(clip);
load(preprocessdata,'ptorder','guidepts','numframes','ptdatmm',...
    'guideptind','arr','numpsplineseg','opt')

% Load process data from file if possible
analysisdatapath = [dropboxpath 'octopus-3d-tracking' filesep 'temp_MATLAB' filesep 'segment_analysis_data_' char(clip) '_' opt.method '.mat'];

if exist(analysisdatapath,"file")==0

    %% Calculate properties of each segment

    % Initialize arrays
    numpts = numel(ptorder)-numel(guidepts);       % Number of tracked (key)points
    numseg = numpts-1;     % Number of segments
    seglenlin = nan(numframes,numseg);  % Linear segment length
    seglenspline = seglenlin;

    ptdatmmNG = ptdatmm(~guideptind,:,:);

    % Loop through frames
    for i=1:numframes
        % Go through segments
        for j=1:numseg
            % Linear segment length
            seglenlin(i,j) = sqrt(sum(diff(ptdatmmNG(j:j+1,:,i),1,1).^2));
        end
    end

    % Define time axis
    timax = (0:numframes-1)/framerate;
    if exist('armlift','var')
        timax = ((0:numframes-1)-(armlift-1))/framerate;
    end

    %% Calculate cumulative distance along points on arm & reinterpolate
    % Make all points equidistant, keeping key points stationary

    splpts = arr;
    warning('Calculating cumulative distance relative to tracked points')

    % Calculate average cumdist through clip
    avcumdist = [0; cumsum(median(sqrt(sum(diff(splpts,1,1).^2,2)),3,'omitnan'))];
    cumdistorig = avcumdist;
    % cumdist = linspace(0,max(avcumdist),interpnum);
    cumdist = 0:min(diff(avcumdist)):max(avcumdist);
    interpnum = numel(cumdist);

    % Fake old spacing to keep keypoints in place
    oldcumdistarr = repmat(cumdist,size(arr,3),1);
    % Interpolate
    newsplpts = nan(interpnum,3,numframes);
    for i=1:numframes
        if sum(sum(~isnan(splpts(:,:,i)))<2>0)
            % Don't worry, just don't fill in any new values if there aren't enough
            % spline points
        else
            temp = ~(sum(isnan(splpts(:,:,i)),2)>0);   % Rows with non-NaN values
            newsplpts(:,:,i) = interp1(avcumdist(temp),splpts(temp,:,i),cumdist,'pchip',NaN);
        end
    end
    splpts = newsplpts;


    %% Calculate angles and stretch

    numsplpts = size(splpts,1);
    spltotang = zeros(numframes,numsplpts-1);

    % splpts = arrGeo;

    xzvec = [0,1,0];    % Normal vector to horizontal plane
    % Loop through frames
    for i=1:numframes
        % Get difference between point locs
        dd = diff(splpts(:,:,i),1);
        % Define plane from vertical axis and previous vector
        for k=2:numsplpts-1
            %     refplane = cross(dd(k-1,:),xzvec);  % Normal vector to vertical plane through previous vector
            % Get relative 'yaw'
            %     splyaw(i,k) = asin(sum(refplane.*dd(k,:),2)./(sqrt(sum(refplane.^2))*sqrt(sum(dd(k,:).^2,2))))';
            % Total angle
            dotprod = sum(dd(k-1,:).*dd(k,:))/(sqrt(sum(dd(k-1,:).^2))*sqrt(sum(dd(k,:).^2)));
            %     if abs(dotprod-1)<1e-7  % Avoid generating complex numbers for aligned vectors
            %         spltotang(i,k) = 0;
            %     else
            %         spltotang(i,k) = acos(dotprod);
            %     end
            if abs(dotprod-1)<1e-7  % Avoid generating complex numbers for aligned vectors
                spltotang(i,k) = 0;
            else
                [R,~] = getSphereFrom3Points(splpts(k-1,:,i),splpts(k,:,i),splpts(k+1,:,i));
                spltotang(i,k) = 1/R;
            end
        end
    end



    %% Plot angles along entire arm; total curvature
    %
    ptinds = [1:(size(arr,1)-1)/(numel(ptorder)-numel(guidepts)-1):size(arr,1)];


    %% Plot curvature surface

    binframes = 10;        % 1 = no binning
    binarm = 0.02;      % Fraction of arm for bin size for curvature

    if ~exist('curvframestart',"var")
        curvframestart = 1;
    end
    if ~exist('curvframeend',"var")
        curvframeend = 1e9;
    end

    curvdat = rad2deg(spltotang);         % Define a curvature parameter to use
    % Smooth just a tiny bit (we are calculating bend radius for very small
    % steps after all
    curvdat = movmean(curvdat,5,2,'omitnan');
    curvframeend = min(curvframeend,size(curvdat,1));

    % Do trimming and binning
    curvdattemp = curvdat(curvframestart:curvframeend,:);
    curvdatbinned = curvdattemp;
    binnum = floor(((1:size(curvdattemp,1))-1)/binframes)+1;
    abinnum = floor((1:size(curvdattemp,2))/(size(curvdattemp,2)*binarm))+1;
    for i=1:max(binnum)
        temp = binnum == i;
        tbin = mean(curvdattemp(temp,:),1,"omitnan");   % Bin in time
        % Bin along arm
        for j = 1:max(abinnum)
            temp2 = abinnum == j;
            vbin(1,temp2) = mean(tbin(1,temp2),2,"omitnan");
        end
        curvdatbinned(temp,:) = repmat(vbin,sum(temp),1);
    end

    % Make array that relates tracked point locations to cumulative distance
    % ptdistarr = oldcumdistarr(:,1:numpsplineseg:end);
    ptdistarr = repmat(avcumdist(1:numpsplineseg:end)',numframes,1);


    %% Show stretching with similar graph as curvature surface above

    alldists = permute(sqrt(sum(diff(arr,1,1).^2,2)),[3,1,2]);
    segdists = nan(numframes,numseg);
    for i=1:numseg
        temp1 = (i-1)*numpsplineseg+1;
        segdists(:,i) = sum(alldists(:,temp1:temp1+numpsplineseg-1),2);
    end

    % Do trimming and binning
    segdiststemp = segdists(curvframestart:curvframeend,:);
    segdistsbinned = segdiststemp;
    binnum = floor(((1:size(segdiststemp,1)-1))/binframes)+1;
    for i=1:max(binnum)
        tempnn = binnum == i;
        segdistsbinned(tempnn,:) = repmat(mean(segdiststemp(tempnn,:),1,"omitnan"),sum(tempnn),1);
    end

%     % Normalize
%     segdistsbinned = 100*segdistsbinned./min(segdistsbinned,[],1,'omitnan')-100;

    segdistbig = curvdatbinned*0;
    for i = 1:numseg
        %     xcvarr = round(ptdistarr(1,i))+1:round(ptdistarr(1,i+1));
        %     segdistbig(:,xcvarr) = repmat(segdists(:,i),1,numel(xcvarr));
        xcvarr = cumdist>=ptdistarr(1,i)&cumdist<ptdistarr(1,i+1);
        segdistbig(:,xcvarr) = repmat(segdistsbinned(:,i),1,sum(xcvarr));
%         if clip=="O15_1611_19083_L3"
%             if i==numseg
%                 segdistbig(:,xcvarr) = repmat(segdistsbinned(:,i)*NaN,1,sum(xcvarr));
%                 warning('Removing row of data for this clip')
%             end
%         end
    end

    %% Get location of max curvature along arm 
    curvPeakInd = nan;
    curvPeakHeight = nan;
    % comment out the line below section if determining armtrimloc
    [curvPeakInd, curvPeakHeight] = curvPeaks(curvdat,armtrimloc);

    %% Save results
    save(analysisdatapath,'cumdist',...
        'timax','curvdat','curvdatbinned','ptdistarr','segdistbig',...
        'splpts','curvPeakInd','curvPeakHeight','segdists');

end
end

%% Functions

function ang = angdif(ang1,ang2)
ang = ang2-ang1;
temp11 = abs(ang)>pi;
ang(temp11) = sign(ang(temp11)).*(2*pi-abs(ang(temp11)));
end

function [peakInd, peakHeight] = curvPeaks(curvdat,cutoffind)
% Cutoff ind is the index above which to ignore entries
% cutoffind = 400;
% If cutoff <= 1, it is assumed to be a relative value
if cutoffind<= 1
    cutoffind = round(cutoffind*size(curvdat,2));
end
cutoffind = min(cutoffind,size(curvdat,2));
curvdat(:,cutoffind:end) = [];
peakInd = nan(size(curvdat,1),1);
peakHeight = peakInd;
firstloop = true;

% Loop through rows (time)
for rrr = 1:size(curvdat,1)
    % Select the current curve, apply a median filter
    thiscurve = movmean(curvdat(rrr,:),[2 2]);
    % Find max
    [maxH,ind] = max(thiscurve);


    if maxH>2

        % Crop a small range around max to find center of peak
        windowsize = 20;
        minind = max(1,ind-windowsize/2);
        maxind = min(numel(thiscurve),ind+windowsize/2);
        inds = minind:maxind;
        cropcurve = thiscurve(inds);

        % Fit Gaussian to small section to find peak
        gaussEqn = 'a*exp(-2*(x-b)^2/c^2)+d';
        temp=isnan(cropcurve);
        tempcurve=cropcurve(~temp);
        indsb=inds(~temp);
        startpoints = [max(tempcurve)-min(tempcurve) indsb(round(end/2)) windowsize/2 0];
        f1 = fit(indsb',tempcurve',gaussEqn,'Start',startpoints,...
            'Robust','Bisquare','Display','off');
        gaussfit = f1.a*exp(-2*(indsb-f1.b).^2/f1.c^2)+f1.d;
        peakInd(rrr) = f1.b;
        peakHeight(rrr) = f1.a;
        figure(16)
        xax = (1:numel(thiscurve));
        if firstloop
            plcurve = plot(xax,thiscurve);
            hold on
            plfit = plot(xax(indsb),gaussfit,'-r');
            hold off
            pltitle = title(num2str(rrr));
            firstloop = false;
        else
            plcurve.XData = xax;
            plcurve.YData = thiscurve;
            plfit.XData = xax(indsb);
            plfit.YData = gaussfit;
            pltitle.String = num2str(rrr);
        end
        pause(0.05);
    end

end

end