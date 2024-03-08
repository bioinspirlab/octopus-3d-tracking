% Calculate geodesic distances between points
% function [pts,dist] = calcPath(depthmap,pt1xy,pt2xy,numpsplineseg,opt)

%% Use for testing
clear all
load geotest.mat

mapData = load("dMapCityBlock.mat");
omap = mapData.omap;
omap = occupancyMap;

omap.FreeThreshold = 0.5;
[XX,YY] = meshgrid(1:100,1:100);
ZZ = sqrt((XX-60).^2+(YY-40).^2);

updateOccupancy(omap,[XX(:),YY(:),ZZ(:)],0.6)
inflate(omap,1)
ss = stateSpaceSE3([0 220;0 220;0 100;inf inf;inf inf;inf inf;inf inf]);
sv = validatorOccupancyMap3D(ss, ...
     Map = omap, ...
     ValidationDistance = 0.1);
planner = plannerRRTStar(ss,sv, ...
          MaxConnectionDistance = 50, ...
          MaxIterations = 1000, ...
          GoalReachedFcn = @(~,s,g)(norm(s(1:3)-g(1:3))<1), ...
          GoalBias = 0.1);
start = [40 180 25 0.7 0.2 0 0.1];
goal = [150 33 35 0.3 0 0.1 0.6];

rng(1,"twister");
[pthObj,solnInfo] = plan(planner,start,goal);
show(omap)
axis equal
view([-10 55])
hold on
% Start state
scatter3(start(1,1),start(1,2),start(1,3),"g","filled")
% Goal state
scatter3(goal(1,1),goal(1,2),goal(1,3),"r","filled")
% Path
plot3(pthObj.States(:,1),pthObj.States(:,2),pthObj.States(:,3), ...
      "r-",LineWidth=2)

ererere
depthmap = double(depthmap);

map3D = occupancyMap3D;
% [xGround,yGround,zGround] = meshgrid(0:100,0:100,0);
xGround = depthmap(:,:,1);
yGround = depthmap(:,:,2);
zGround = depthmap(:,:,3);
nanlocs = isnan(zGround(:));
zGround(nanlocs) = 0;
% xyzGround = [xGround(~nanlocs) yGround(~nanlocs) zGround(~nanlocs)];
xyzGround = [xGround(:) yGround(:) zGround(:)];
occval = 0;
setOccupancy(map3D,xyzGround,occval);


obstr = [xGround(nanlocs) yGround(nanlocs) zGround(nanlocs)];
updateOccupancy(map3D,obstr,1)
map3D.FreeThreshold = 0.5;
show(map3D)

ss = stateSpaceSE3([min(xGround(:)) max(xGround(:));min(yGround(:)) max(yGround(:));min(zGround(:)) max(zGround(:));inf inf;inf inf;inf inf;inf inf]);
sv = validatorOccupancyMap3D(ss, ...
     Map = map3D, ...
     ValidationDistance = 0.1);
planner = plannerRRTStar(ss,sv, ...
          MaxConnectionDistance = 50, ...
          MaxIterations = 1000, ...
          GoalReachedFcn = @(~,s,g)(norm(s(1:3)-g(1:3))<1), ...
          GoalBias = 0.1);

sz = size(depthmap);
[XX,YY] = meshgrid(1:sz(2),1:sz(1));
% start = [interp1(1:sz(1),depthmap(1,:,1),pt1xy(1)) interp1(1:sz(1),depthmap(:,1,2),pt1xy(2)) interp2(XX,YY,depthmap(:,:,3),pt1xy(1),pt1xy(2)) 0.2 0 0.1 0];
% goal = [interp1(1:sz(1),depthmap(1,:,2),pt2xy(1)) interp1(1:sz(1),depthmap(:,1,2),pt2xy(2)) interp2(XX,YY,depthmap(:,:,3),pt2xy(1),pt2xy(2)) 0.2 0 0.1 0];
start = [permute(depthmap(round(pt1xy(1)),round(pt1xy(2)),:),[2 3 1]) 0.2 0 0.1 0];
goal = [permute(depthmap(round(pt2xy(1)),round(pt2xy(2)),:),[2 3 1]) 0.2 0 0.1 0];
[pthObj,solnInfo] = plan(planner,start,goal);

show(map3D)
axis equal
view([-10 55])
hold on
% Start state
scatter3(start(1,1),start(1,2),start(1,3),"g","filled")
% Goal state
scatter3(goal(1,1),goal(1,2),goal(1,3),"r","filled")
% Path
plot3(pthObj.States(:,1),pthObj.States(:,2),pthObj.States(:,3), ...
      "r-",LineWidth=2)
% plannermap = isnan(depthmap(:,:,3));

% end