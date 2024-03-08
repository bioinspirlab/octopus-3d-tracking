% Calculate geodesic distances between points
function [pts,dist] = calcGeodesic(depthmap,pt1xy,pt2xy,numpsplineseg,opt)

% Use for testing
% load geotest.mat

if nargin<5
    opt.resize = true;  % Default: resize
    opt.resizefactor = 'auto';
end

%% Adaptive resize for speed
if opt.resize
    if ischar(opt.resizefactor)
        % Automatically determine suitable resize factor
        eucliddist = sqrt(sum((pt1xy-pt2xy).^2));
        rsfact = min(3*numpsplineseg/eucliddist,1);
    else
        rsfact = opt.resizefactor;
    end
    depthmap = imresize(depthmap,rsfact);
    pt1xy = (pt1xy+1)*rsfact;
    pt2xy = (pt2xy+1)*rsfact;
end

%%
figure(8)
% depthmap(:,:,3) = medfilt2(depthmap(:,:,3),[15,15]);
% depthmap(:,:,4) = [];
[ysize,xsize,~] = size(depthmap);
D = nan(ysize,xsize); % Distance map
Dsrc = D;   % Map of source pixels to track path
reached = false;
[X,Y] = meshgrid((1:size(depthmap,1)),(1:size(depthmap,2)));
X1 = X-1-pt1xy(1);
Y1 = Y-1-pt1xy(2);

R = sqrt(X1.^2+Y1.^2);
eucliddist = sqrt(sum((pt1xy-pt2xy).^2));
% Cut off R outside twice this distance
R(R>4*eucliddist) = NaN;

% Make an order list of entries to calculate
[~,ind] = min(R(:));
updatelist = [ind, nan];

% Interpolant
Fx = griddedInterpolant(X',Y',depthmap(:,:,1)');
Fy = griddedInterpolant(X',Y',depthmap(:,:,2)');
Fz = griddedInterpolant(X',Y',depthmap(:,:,3)');

refX = Fx(pt1xy);
refY = Fy(pt1xy);
refZ = Fz(pt1xy);

eucliddistmm = sqrt(sum(([refX,refY,refZ]-[Fx(pt2xy),Fy(pt2xy),Fz(pt2xy)]).^2));

i=1;
Dsz = size(D);
while size(updatelist,1)>0
    ind = updatelist(1,1);
    xx = X(ind);
    yy = Y(ind);
    curr = depthmap(yy,xx,:);
    if sum(isnan(curr))>0
        disp('ff')
    end

    if i==1
        if sum(isnan(curr))>0
            pts = [nan, nan, 0];
            warning('No valid data on first point in geodesic analysis.')
            break
        end
        D(yy,xx) = sqrt(sum((permute([refX,refY,refZ],[3,1,2])-curr).^2));
    end

%     ylocs = max(yy-1,1):min(yy+1,ysize);
%     xlocs = max(xx-1,1):min(xx+1,xsize);
     ylocs = max(yy-2,1):min(yy+2,ysize);
    xlocs = max(xx-2,1):min(xx+2,xsize);
    temp = sqrt(sum((depthmap(ylocs,xlocs,:)-curr).^2,3))+D(yy,xx); % Distances

    %     temp2=temp>=D(ylocs,xlocs);  % Values we do NOT want to paste
    temp2=temp>=D(ylocs,xlocs)|isnan(temp); % Values we do NOT want to paste

    if sum(temp2,[1 2])>0
        % If we're picky
        [c,v]=meshgrid(xlocs,ylocs);
        linearidx = sub2ind(Dsz,v(~temp2),c(~temp2));
        D(linearidx) = temp(~temp2);
        Dsrc(linearidx) = ind;
        if reached
            aa = temp(~temp2);
            updatelist = [updatelist; [linearidx(aa<destlev),aa(aa<destlev)]];
        else
            updatelist = [updatelist; [linearidx,temp(~temp2)]];
        end
    else
        % If we want to write all entries
        D(ylocs,xlocs) = temp;
        Dsrc(ylocs,xlocs) = ind;    % Tell replaced pictures their connector
    end
    updatelist(updatelist(:,1)==ind,:)=[];

    % Re-sort due to new entries; small distances first
    if mod(i,500)==0
        [~,temp]=sort(updatelist(:,2));
        updatelist = updatelist(temp,:);
    end

    if mod(i,4000)==1
        % Clean up list items getting out of defined bounds
        updatelist(isnan(R(updatelist(:,1))),:)=[];

        % Determine if we've reached our destination
        destlev = D(round(pt2xy(2)),round(pt2xy(1)));
        if ~isnan(destlev)
            reached = true;
            % Clear hopeless points
            updatelist(updatelist(:,2)>destlev,:) = [];
        end

        if i==1
            figure(8)
            imm = imagesc(D);
            axx = gca;
        else
            imm.CData = D;
        end

        if i>10000000
            warning('Geodesic calculations timed out');
            break
        end
        drawnow;
        pause(0.002);

        if i>numel(depthmap)
            break
        end
    end
    i=i+1;
end

if i==1
    warning('Terminating geodesic analysis: could not start')
    pts = [pt1xy;pt2xy];
    dist = [0; eucliddistmm];

    return

end

imm.CData = D;
hold(axx,"on");

% Now determine path, going backwards

% Find first and last point first
try
pts = Dsrc(round(pt2xy(2)),round(pt2xy(1)));
if isnan(pts)
    for k = 1:10 % Look in increasing box nearby for valid points
        if isnan(pts)
    temp=D(round(pt2xy(2))+(-k:k),round(pt2xy(1))+(-k:k));
    [~,temp]=min(temp,[],'all');
    [r,c] = ind2sub([2*k+1 2*k+1],temp);
    pts = Dsrc(round(pt2xy(2))-k-1+r,round(pt2xy(1))-k-1+c);
        end
    end
end
for i=1:500
    pts = [pts Dsrc(pts(end))];
    if R(pts(end))<1
        break
    end
end
catch
    warning('Could not find valid path; terminating geodesic analysis')
    pts = [pt1xy;pt2xy];
    dist = [0; eucliddistmm];

    return
end
[row,col] = ind2sub(Dsz,pts);
pts = [col',row'];
% pts = flipud(pts);
% pts(:,3)=[0; cumsum(sqrt(sum(diff([Fx(pts(:,1:2)),Fy(pts(:,1:2)),Fz(pts(:,1:2))],1,1).^2,2)))];
% plot(axx,pts(:,1),pts(:,2),'-b');
% 
% % 
% % 
% figure(8)
% hold off
%             imm = imagesc(D);
%             axx = gca;
%             hold on
% pts = pt2xy;
% [FX,FY] = gradient(D);
% ang = griddedInterpolant(X',Y',atan2(FY,FX));
% step = 1;
% i=1;
% while true
%     
%     xx=pts(end,1);
%     yy=pts(end,2);
% 
%     % Terminate when close to pt1
%     if sqrt(sum((pts(end,1:2)-pt1xy).^2))<3
%         break
%     end
%     % Otherwise, add point to list
%     currang = ang(yy,xx);
%     if isnan(currang)
%         if size(pts,1)==1
%             currang = ang(round(pt2xy(2)),round(pt2xy(1)));
%             pts(end+1,1:2) = [xx-step*cos(currang),yy-step*sin(currang)];
%         else
%             currang = ang(pts(end-1,2),pts(end-1,1));
%             tempstep = 1:20;
%             temp = ang(yy-tempstep*sin(currang),xx-tempstep*cos(currang));
%             ind = find(~isnan(temp),1,'first');
%             if ~isempty(ind)
%                 pts(end+1,1:2) = [xx-tempstep(ind)*cos(currang),yy-tempstep(ind)*sin(currang)];
%             else
%                 % Give up
%                 warning('Terminating geodesic analysis due to a problem')
%                 pts = [pt1xy;pt2xy];
%                 dist = [0; eucliddistmm];
% 
%                 return
%             end
%         end
%     else
% 
%         pts(end+1,1:2) = [xx-step*cos(currang),yy-step*sin(currang)];
%         if size(pts,1)>1e5
%             warning(['Terminating geodesic analysis due to a problem:' ...
%                 ' trapped in local minimum'])
%             pts = [pt1xy;pt2xy];
%             dist = [0; eucliddistmm];
% 
%             return
%         end
%     end
%     hold on
%     plot(pts(end,1),pts(end,2),'.r')
%     pause(1)
%     i=i+1;
% 
% end
% pts(end+1,:) = pt1xy;

pts = flipud(pts);
pts(:,3)=[0; cumsum(sqrt(sum(diff([Fx(pts(:,1:2)),Fy(pts(:,1:2)),Fz(pts(:,1:2))],1,1).^2,2)))];
plot(axx,pts(:,1),pts(:,2),'-r');
hold(axx,"off");

% If there is a large jump in point distance (likely due to z)
temp = abs(diff(pts(:,3)));
if sum(temp>max(50*median(temp),(median(temp)+10*std(temp))))>0
    warning('Large jump detected; reverting to euclidean');
    pts = [pt2xy eucliddistmm; nan nan eucliddistmm/2];
end

% Output points; reduce number to make smoother
% ptlocs = round(linspace(1,size(pts,1),max(size(pts,1)/10,min(3,size(pts,1)))));
% pts = pts(ptlocs,:);
interploc = (0:numpsplineseg-1)*pts(end,3)/numpsplineseg;

pts = interp1(pts(:,3),pts,interploc);
dist = pts(:,3);
pts(:,3) = [];
if opt.resize
    pts = pts/rsfact-1;
end

% disp(['Geodesic distance = ' num2str(max(dist))])
% disp(['Euclidean distance = ' num2str(eucliddistmm)])
euclidratio = max(dist)/eucliddistmm;
if euclidratio>5
    disp(['Euclidean distance ratio = ' num2str(euclidratio)])
end

end

