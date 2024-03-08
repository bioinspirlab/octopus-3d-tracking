% Calculate spline points
function [splinex,spliney,splinezmm] = surfaceSpline(ptdat,ii,ptorder,numpsplineseg,zmm,opt)
if nargin<6
    opt.method = 'makima';
elseif ~isfield(opt,'method')
    opt.method = 'makima';
end
if isfield(opt,'guidepts')
    guidepts = opt.guidepts;
end

validdat = permute(ptdat(ii,:,ptorder),[3 2 1]);
validdat(:,3) = zmm(ptorder);
badpts = isnan(validdat(:,1));
% badptlocs=find(badpts);
goodptlocs = find(~badpts);
validdat(badpts,:) = [];

% Ignore guide points in spline lengths
keyptorder = ptorder(~ismember(ptorder,guidepts));
keypointdat = find(~ismember(ptorder(goodptlocs),guidepts));

numsplpts = numpsplineseg*(numel(ptorder)-numel(guidepts)-1)+1; % Number of spline points

if numel(goodptlocs)<2 || numel(keypointdat)<2
    % Need at least two points for interpolation
    splinex = nan(numsplpts,1);
    spliney = splinex;
    splinezmm = splinex;
    if numel(goodptlocs)==1
        ptloc = find(keyptorder==ptorder(goodptlocs(end)));
        splinex((ptloc-1)*numpsplineseg+1) = validdat(1,1);
        spliney((ptloc-1)*numpsplineseg+1) = validdat(1,2);
        splinezmm((ptloc-1)*numpsplineseg+1) = validdat(1,3);
    end
    return
end

% % Set up interpolation axis - old way without guide points
% CS = cat(1,0,cumsum(sqrt(sum(diff(validdat(:,1:2),[],1).^2,2))));
% CSint = nan(numpsplineseg*(numel(ptorder)-1)+1,1);
% zval = CSint;
% for i=1:size(CS,1)-1
%     jumpsz = goodptlocs(i+1)-goodptlocs(i);
%     intv = (CS(i+1)-CS(i))/(numpsplineseg*jumpsz);
%     CSint((goodptlocs(i)-1)*numpsplineseg+(1:numpsplineseg*jumpsz)) = ...
%         (0:jumpsz*numpsplineseg-1)*intv+CS(i);
%     if ~isnan(validdat(i,3)+validdat(i+1,3))
%         zval((goodptlocs(i)-1)*numpsplineseg+(1:numpsplineseg*jumpsz)) = 1;
% 
%     end
% end
% CSint((goodptlocs(end)-1)*numpsplineseg+1) = CS(end);
% if ~isnan(validdat(end,3))
%     zval((goodptlocs(end)-1)*numpsplineseg+1) = 1;
% end

% Set up interpolation axis - new way with guide points
% Approach: determine spacing from track points, ignoring guide points
% (they will still be used in interpolation below)
CS = cat(1,0,cumsum(sqrt(sum(diff(validdat(:,1:2),[],1).^2,2))));
CSint = nan(numsplpts,1);
zval = CSint;
CSkey = CS(keypointdat);
for i=1:size(CSkey,1)-1
    currpt = ptorder(goodptlocs(keypointdat(i)));
    numpts = goodptlocs(keypointdat(i+1))-goodptlocs(keypointdat(i));
    numguidepts = sum(ismember(ptorder(goodptlocs(keypointdat(i)):goodptlocs(keypointdat(i+1))),guidepts));
    jumpsz = numpts-numguidepts;
    intv = (CSkey(i+1)-CSkey(i))/(numpsplineseg*jumpsz);
%     CSint((find(keyptorder==currpt)-1)*numpsplineseg+(1:numpsplineseg*jumpsz)) = ...
%         (0:jumpsz*numpsplineseg-1)*intv+CSkey(i);
    CSint((find(keyptorder==currpt)-1)*numpsplineseg+(1:numpsplineseg*jumpsz+1)) = ...
        (0:jumpsz*numpsplineseg)*intv+CSkey(i);
    if ~isnan(validdat(keypointdat(i),3)+validdat(keypointdat(i+1),3))
        zval((find(keyptorder==currpt)-1)*numpsplineseg+(1:numpsplineseg*jumpsz+1)) = 1;
    end
end
CSint((find(keyptorder==ptorder(goodptlocs(keypointdat(end))))-1)*numpsplineseg+1) = CSkey(end);
if ~isnan(validdat(end,3))
    zval((find(keyptorder==ptorder(goodptlocs(keypointdat(end))))-1)*numpsplineseg+1) = 1;
end

% Do the interpolation
znan = ~isnan(validdat(:,3));
% lastz = find(znan,1,'last');
dd = interp1(CS, validdat(:,1:2), CSint,opt.method);
splinex = dd(:,1);
spliney = dd(:,2);
% splinezmm = dd(:,3);
if numel(validdat(znan,3))>1
    splinezmm = interp1(CS(znan), validdat(znan,3), CSint.*zval,opt.method);
else
    splinezmm = nan(numsplpts,1);
    splinezmm(~isnan(zval)) = validdat(znan,3);
end
end