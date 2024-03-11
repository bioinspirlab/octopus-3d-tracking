%% Information

% EyeRIS octopus analysis from DLTdv points
% Converts points to XYZ points in mm using the depthmap data, performs
% spline interpolation between 3D points, and exports a preview clip

function processdatapath = octo_PreProcess(clip)

%% Configure
matpath = octo_InitializeData(clip);
load(matpath);

maxgapfill = 100;
inspect = false;
opt.method = 'spline';

%% Load existing file if available

% Load process data from file if possible
processdatapath = [dropboxpath 'octopus-3d-tracking' filesep 'temp_MATLAB' filesep 'preprocessed_data_' char(clip) '_' opt.method '.mat'];
if exist(processdatapath,"file")==0

    %% Load point data
    rawdata = importdata(datapath);
    numpoints = size(rawdata.data,2)/2;
    if ~exist("numframes","var")
        numframes = size(rawdata.data,1);
    end

    ptdat = zeros(numframes,2,numpoints);
    ptdatmm = nan(numpoints,3,numframes);
    for i=1:numpoints
        ptdat(:,1,i) = rawdata.data(1:numframes,(i-1)*2+1); % x data
        ptdat(:,2,i) = rawdata.data(1:numframes,(i-1)*2+2); % y data
    end

    % Correct a weird time offset issue with the 7318 data
    switch clip
        case {"octo_1611_7318_L1","octo_1611_7318_L2"}
            warning('Special data correction happens for this clip');
            ptnumcorr = [2,3,4,5,6,10:20];
            ptdat(500:end-1,:,:) = ptdat(501:end,:,:);
    end

    if exist("startnum","var")==0
        temp = strsplit(clip,'_');
        startnum = str2double(temp(end));
    end

    keyptorder = ptorder;
    for i=1:numel(guidepts)
        keyptorder(keyptorder==guidepts(i)) = [];
    end

    %% Loop through data, getting Z position from depthmap

    firstloop = true;
    reverseStr = [];
    for i=1:numframes
        newstr = ['Preprocessing... frame ', num2str(i), '/', num2str(numframes)];
        fprintf([reverseStr, newstr]);
        reverseStr = repmat(sprintf('\b'), 1, length(newstr));
        try % Allow continuation even if one loop fails

            % Only continue if there is data for this frame
            if sum(~isnan(sum(ptdat(i,:,keyptorder),2)))>0

                imname_depth = [orthopath,num2str(startnum-1+i,'%04d'),'.tif'];
                temp_depth = imread(imname_depth);

                % Initialize arrays on first loop
                if firstloop
                    % Horizontal axis
                    temp = double(median(temp_depth(:,:,3),1,'omitnan'));
                    xpix = 1:size(temp_depth,2);
                    ind = ~isnan(temp);
                    xarr = interp1(xpix(ind),temp(ind),xpix,'linear','extrap');

                    % Vertical axis
                    temp = double(median(temp_depth(:,:,2),2,'omitnan'));
                    ypix = 1:size(temp_depth,1);
                    ind = ~isnan(temp);
                    yarr = interp1(ypix(ind),temp(ind),ypix,'linear','extrap');

                    % Meshgrid for interpolation on depthmap
                    [X,Y] = meshgrid(xpix,ypix);

                    firstloop = false;
                end

                % Generate depthmap
                [depthmap,opt] = octo_ProcessDepthMap(i,orthopath,startnum,opt);

                xloc = ptdat(i,1,:);
                yloc = ptdat(i,2,:);
                ptdatmm(:,1,i) = permute(interp1(xpix,xarr,xloc),[3,1,2]);
                ptdatmm(:,2,i) = permute(interp1(ypix,yarr,yloc),[3,1,2]);
                ptdatmm(:,3,i) = permute(interp2(X,Y,depthmap,xloc,yloc),[3,1,2]);
            end

        catch ME
            warning(['Loop ' num2str(i) ' failed. ' ME.message])
            reverseStr = [];
        end
    end

    %% Data cleanup
    % Trim data depending on clip, removing points in a specific window of
    % frames

    for i=1:size(erasethese,1)
        ptnum = erasethese(i,1);
        erasestart = erasethese(i,2);
        eraseend = erasethese(i,3);
        ptdatmm(ptnum,3,erasestart:eraseend) = NaN;
        %     ptdat(erasestart:eraseend,:,ptnum) = NaN;
    end

    if startfr>1 || endfr<size(ptdatmm,3)
        ptdat = ptdat(startfr:endfr,:,:);
        ptdatmm = ptdatmm(:,:,startfr:endfr);
    end

    numframes = size(ptdatmm,3);
    
    ptorderNG = ptorder(~ismember(ptorder,guidepts));  % Point numbers in order, excluding guide pts


    %% Median filter to correct outliers in reference point z data

    ptdatmmnans = isnan(ptdatmm);
    ptdatmm(refpts,:,:) = movmedian(ptdatmm(refpts,:,:),5,3,'omitnan');
    ptdatmm(ptdatmmnans) = NaN;


    %% Visual check of reference data (if requested)
    if inspect
        figure(98)
        refdat = ptdatmm(refpts,:,:);
        refnum = cellstr(num2str(refpts'));
        plot(squeeze(refdat(:,1,:))')
        hold off
        legend(refnum)
        title('Reference points; X coordinate')
        plot(squeeze(refdat(:,2,:))')       % MARK HERE
        legend(refnum)
        title('Reference points; Y coordinate')
        plot(squeeze(refdat(:,3,:))')       % MARK HERE
        legend(refnum)
        title('Reference points; Z coordinate')
        clear refdat    %       MARK HERE
    end

    %% Data prep: subtract background motion and remove reference points

    % Subtract baseline movement from pixel location array
    bgmove_px = median(diff(ptdat(:,:,refpts),1,1),3,'omitnan');
    bgmove_px(isnan(bgmove_px)) = 0;
    bgmove_px = cumsum(bgmove_px,1);
    ptdat(2:end,:,:) = ptdat(2:end,:,:)-bgmove_px;

    % Subtract baseline movement from 3D (mm) arrays
    intx = median(diff(double(xarr)));
    inty = median(diff(double(yarr)));
    bgmove_mm = permute(bgmove_px.*[intx inty],[3,2,1]);
    basedifz = median(diff(ptdatmm(refpts,3,:),1,3),1,'omitnan');
    basedifz(isnan(basedifz)) = 0;
    basedifz = cumsum(basedifz,3);
    basedifz = movmedian(basedifz,7,3,'omitnan'); % Minor smoothing of depth data
    bgmove_mm(:,3,:) = basedifz;
    ptdatmm(:,:,2:end) = ptdatmm(:,:,2:end)-bgmove_mm;

    % Put points in order, remove background reference points
    ptdat = ptdat(:,:,ptorder);
    ptdatmm = ptdatmm(ptorder,:,:);
    ptnames = cellstr(num2str(ptorder'));
    guideptind = ismember(ptorder,guidepts);
    guidepts = find(guideptind);
    ptorderorig = ptorder;
    ptorder = 1:numel(ptorder);


    %% Smooth data (remove outliers, smooth data, fill small gaps)

    if inspect
        figure(98)
        ppp = plot(squeeze(ptdatmm(:,3,:))'+20*(1:numel(ptorder)));
        legend(ptnames)
        hold on
    end
    % Remove outliers in depth data
    if ~exist('thresh','var')
        thresh = 1;     % Multiple of standard deviation to consider outlier
    end
    zref = movmedian(ptdatmm(:,:,:),50,3,'omitnan');
    dev = abs(zref-ptdatmm(:,:,:));
    devstd = std(dev,0,3,'omitnan');
    devstd(devstd<5)=5;
    inval = dev>thresh*devstd;    % Invalid indices
    inval(:,1:2,:) = false;
    ptdatmm(inval) = NaN;
    if inspect
        figure(98)
        clrs = {ppp.Color};
        ppp2 = plot(squeeze(ptdatmm(:,3,:))'+20*(1:numel(ptorder)),'linewidth',4);
        for i=1:numel(ppp)
            ppp2(i).Color = clrs{i};
        end
        hold off
        legend([ptnames;ptnames])
        title('Removing outliers (z position, each offset by 20 mm for clarity)')
        pause(0.1)      % MARK HERE
    end

    % Smooth data
    xysmooth = 20;   % Smoothing window in image plane
    zsmooth = 60;   % Smoothing window of depth values
    ptdatnans = isnan(ptdat);
    ptdat(:,1,:) = movmean(ptdat(:,1,:),xysmooth,1,'omitnan');
    ptdat(:,2,:) = movmean(ptdat(:,2,:),xysmooth,1,'omitnan');
    ptdat(ptdatnans) = nan;
    % ptdatmmnans = isnan(ptdatmm(:,1,:));
    % ptdatmmnans = repmat(ptdatmmnans,1,3,1);
    ptdatmmnans = isnan(ptdatmm);
    ptdatmm(:,1,:) = movmean(ptdatmm(:,1,:),xysmooth,3,'omitnan');
    ptdatmm(:,2,:) = movmean(ptdatmm(:,2,:),xysmooth,3,'omitnan');
    ptdatmm(:,3,:) = movmean(ptdatmm(:,3,:),zsmooth,3,'omitnan');
    ptdatmm(ptdatmmnans) = nan;

    if inspect
        figure(98)
        hold on
        ppp3 = plot(squeeze(ptdatmm(:,3,:))'+20*(1:numel(ptorder)),'-k','LineWidth',1);
        hold off
        pause(0.1)      % MARK HERE
    end

    % Fill gaps in point data smaller than maxgapfill
    frarray = 1:numframes;
    % Loop through points
    keypts = 1:size(ptdat,3);
    keypts = keypts(~guideptind);
    for pt = keypts     % Loop through key points only, don't interpolate guide points
        % Fix small gaps in XY
        missing = findSmallGaps(sum(ptdat(:,:,pt),2),maxgapfill);
        if sum(missing)>0
            if inspect
                figure(98)
                plot(ptdat(:,:,pt),'linewidth',5); hold on
            end
            ptdat(missing,:,pt) = interp1(frarray(~missing),ptdat(~missing,:,pt),frarray(missing),'pchip');
            ptdatmm(pt,1,missing) = interp1(frarray(~missing),squeeze(ptdatmm(pt,1,~missing)),frarray(missing),'pchip');
            ptdatmm(pt,2,missing) = interp1(frarray(~missing),squeeze(ptdatmm(pt,2,~missing)),frarray(missing),'pchip');
            if inspect
                plot(ptdat(:,:,pt),'linewidth',1);
                title(num2str(ptorderorig(pt)));
                hold off    % MARK HERE
            end
        end
        %     pause(1);

        % Fix small gaps in Z
        missing = findSmallGaps(ptdatmm(pt,3,:),maxgapfill);
        if sum(missing)>0
            if inspect
                plot(squeeze(ptdatmm(pt,3,:)),'linewidth',5); hold on
            end
            ptdatmm(pt,3,missing) = interp1(frarray,squeeze(ptdatmm(pt,3,:)),frarray(missing),'pchip');
            if inspect
                plot(squeeze(ptdatmm(pt,3,:)),'linewidth',1);
                title(['Z coordinate, pt. ' num2str(ptorderorig(pt))]);
                hold off    % MARK HERE
            end
        end

    end

    %% Calculate spline with smoothed data

    numframes = size(ptdatmm,3);
    intx = median(diff(double(xarr)));
    inty = median(diff(double(yarr)));
    opt.guidepts = guidepts;
    % Spline interpolation
    numpsplineseg = 30;
    numpspline = numpsplineseg*(numel(ptorder)-numel(guidepts)-1)+1;
    arr = nan(numpspline,3,numframes);
    xarr = double(xarr);
    yarr = double(yarr);

    for i=1:numframes
        [splinex,spliney,splinezmm] = surfaceSpline(ptdat,i,ptorder,numpsplineseg,ptdatmm(:,3,i),opt);
        splinexmm = (splinex-1)*intx+xarr(1);
        splineymm = (spliney-1)*inty+yarr(1);

        splmm = [splinexmm,splineymm,splinezmm];
        % If z is missing, might as well clear x and y
        znans = isnan(splinezmm);

        % Determine new distances for equal spacing between keypoints
        cumdistall = nan(numel(splinexmm),1);
        cumdistall(~znans) = [0; cumsum(sqrt(sum((diff(splmm(~znans,:),1,1)).^2,2)))];
        keyptloc = 1:numpsplineseg:numel(splinexmm);
        newarr = interp1(keyptloc,cumdistall(keyptloc),1:numel(splinexmm),'linear')';
        if sum(~znans)>1        % If at least one valid entry (there is nothing to do otherwise)
            splmm = interp1(cumdistall(~znans),splmm(~znans,:),newarr,'spline');
        end

        arr(:,:,i) = splmm;
    end

    % Smooth spline (7 frame moving mean)
    arrnans = isnan(arr);
    arr = movmean(arr,7,3,'omitnan');
    arr(arrnans) = nan;

    % Save data
%     save(processdatapath,'ptdat','ptdatmm','arr','clip','framerate','xarr','yarr','startnum')
    save(processdatapath)

end
end

function missing = findSmallGaps(arr,maxgap)
% All missing data in XY (in ptdat, which is used in spline recalc in the next section)
arr = arr(:);   % Work with column array
missing = isnan(arr);

% Do not extrapolate (do not fill nans at beginning or end)
missing(1:find(missing==0,1,"first")) = 0;
missing(find(missing==0,1,"last"):end) = 0;

% Select only gaps smaller than given size
% Inspiration: https://www.mathworks.com/matlabcentral/answers/382011-how-to-count-the-number-of-consecutive-identical-elements-in-both-the-directions-in-a-binary-vecto
d = [true, diff(missing') ~= 0, true];  % TRUE if values change
n = diff(find(d));               % Number of repetitions
Y = repelem(n, n);
missing(Y>maxgap) = false;
end

%% OLD geodesic code
     %             if mode == "geodesic"
        %                 % Complicated and slow geodesic interpolation
        %                 validdat = ptdat(i,:,ptorder);
        %                 badpts = permute(isnan(validdat(1,1,:)),[3,1,2]);
        %                 goodpts = find(~badpts);
        %                 validdat(:,:,permute(isnan(validdat(1,1,:)),[3,1,2])) = [];
        %                 allpts = nan((numel(ptorder)-1)*numpsplineseg+1,2);
        %                 totdist = nan((numel(ptorder)-1)*numpsplineseg+1,1);
        %                 for k = 1:numel(goodpts)-1%numel(ptorder)-1
        %                     %         [pts,dist] = calcGeodesic(temp_depth,ptdat(i,1:2,ptorder(k)),ptdat(i,1:2,ptorder(k+1)));
        %                     %
        %                     temp_depth = xM;
        %                     temp_depth(:,:,2) = yM;
        %                     temp_depth(:,:,3) = depthmap;
        %                     [pts,dist] = calcGeodesic(temp_depth,validdat(:,1:2,k),validdat(:,1:2,k+1),numpsplineseg);
        %
        %                     %                   % Do again, in more detail
        %                     %                     ptind = sub2ind(size(temp_depth),round(pts(:,2)),round(pts(:,1)));
        %                     %                     bw = false(size(temp_depth,1),size(temp_depth,2));
        %                     %                     bw(ptind) = true;
        %                     %                     bw = imdilate(bw,strel('disk',40));
        %                     %                     geo_opt.resize = true;
        %                     %                     geo_opt.resizefactor = 'auto';
        %                     %                     bw = repmat(bw,[1 1 3]);
        %                     %                     temp_depth(~bw) = NaN;
        %                     % [pts,dist] = calcGeodesic(temp_depth,validdat(:,1:2,k),validdat(:,1:2,k+1),numpsplineseg,geo_opt);
        %
        %
        %                     stloc = ((goodpts(k)-1)*numpsplineseg+1);
        %                     allpts(stloc:stloc+size(pts,1)-1,:) = pts;
        %                     totdist(stloc:stloc+size(pts,1)-1,:) = dist;
        %
        %                     %
        %                     %                     allpts = [allpts; pts];
        %                     %                     if k==1
        %                     %                         totdist = dist;
        %                     %                     else
        %                     %                         totdist = [totdist; dist+totdist(end)];
        %                     %                     end
        %                 end
        %                 % Remove duplicates
        %                 %                 allpts(isnan(totdist),:) = [];
        %                 %                 totdist(isnan(totdist)) = [];
        %                 %                 [~,temp] = unique(totdist);
        %                 %                 totdist = totdist(temp);
        %                 %                 allpts = allpts(temp,:);
        %                 %                 geo = interp1(totdist, allpts, linspace(0,totdist(end),numpspline),'makima');
        %                 %                 arrGeo(:,1,i) = interp2(X,Y,xM,geo(:,1),geo(:,2));
        %                 %                 arrGeo(:,2,i) = interp2(X,Y,yM,geo(:,1),geo(:,2));
        %                 %                 arrGeo(:,3,i) = interp2(X,Y,depthmap,geo(:,1),geo(:,2));
        %                 if ~isempty(allpts)
        %                     arrGeo(:,1,i) = interp2(X,Y,xM,allpts(:,1),allpts(:,2));
        %                     arrGeo(:,2,i) = interp2(X,Y,yM,allpts(:,1),allpts(:,2));
        %                     arrGeo(:,3,i) = interp2(X,Y,depthmap,allpts(:,1),allpts(:,2));
        %                 end
        %             end