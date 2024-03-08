function [R,C] = getSphereFrom3Points(p1,p2,p3)
% Calculate sphere radius and center from three points on it
% Inspiration: https://mathworld.wolfram.com/Plane-PlaneIntersection.html

% % Known sphere (from circle in XY plane with radius 3 and center (0,-1,0).
% p1 = [-3, -1, 0];
% p2 = [0, 2, 0];
% p3 = [3/sqrt(2), 3/sqrt(2)-1, 0];

% Ensure we are using row vectors
p1 = p1(:)';
p2 = p2(:)';
p3 = p3(:)';

% Find vectors connecting points (AKA normal vector to bisection planes)
n1 = (p1-p2);     
n2 = (p2-p3);
% Plane reference points
x1 = (p2+p1)/2;
x2 = (p2+p3)/2;

% Check they're not collinear
if sum(cross(n1,n2))==0
    warning('Vectors are collinear; bend radius is infinite')
    R = Inf;
    C = [nan nan nan];
    return
end

% For final vector, use normal vector of plane containing all 3 points.
n3 = cross(n1,n2);
x3 = p3;

C = det([n1;n2;n3])^(-1)*(dot(x1,n1)*cross(n2,n3)+dot(x2,n2)*cross(n3,n1)+dot(x3,n3)*cross(n1,n2));
R = sqrt(sum((C-p1).^2));
% R2 = sqrt(sum((C-[p1;p2;p3]).^2,2))       % Quick test that all distances
% are equal

end