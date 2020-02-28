 function [p1, p2] = max_distance_callipers(points)
%% This function finds the two most distant points on a curve using the rotating callipers algorithm.
%  
%  Input - 2D profile curve (points) with two columns x,y
%  Output - Two points p1, p2 with columns x,y
% 
% This function is dependent on the package geom2d-2018.06.07.
%%  Compute convex hull of the polygon
x = points(:,1);
y = points(:,2);
inds = convhull(x,y);
hull = points(inds,:);
% plot(u(inds),v(inds),'*g',u,v,'.b')

nV = size(hull, 1);                                                         % number of hull vertices
%% Default values
rotated_angle = 0;
max_width = 0;
min_angle = 0;

%Finding points with max and min y-values for an airfoil oriented in the y-direction
[tmp, p_a] = max(hull(:, 2));                                   
[tmp, p_b] = min(hull(:, 2));

caliper_a = [ 1 0];    % Caliper A points along the positive x-axis
caliper_b = [-1 0];    % Caliper B points along the negative x-axis
%% 

while rotated_angle <= pi
    % compute the direction vectors corresponding to each edge
    ind_a2 = mod(p_a, nV) + 1;
    vector_a = hull(ind_a2, :) - hull(p_a, :);
    
    ind_b2 = mod(p_b, nV) + 1;
    vector_b = hull(ind_b2, :) - hull(p_b, :);
    %% 
    % Determine the angle between each caliper and the next adjacent edge in the polygon 
    angle_a = vectorAngle(caliper_a, vector_a);
    angle_b = vectorAngle(caliper_b, vector_b);
    
    % Determine the smallest of these angles
    minAngle = min(angle_a, angle_b);
    
    % Rotate the calipers by the smallest angle
    caliper_a = rotateVector(caliper_a, minAngle);
    caliper_b = rotateVector(caliper_b, minAngle);
    
    rotated_angle = rotated_angle + minAngle;
    %%  compute current width, and update opposite vertex
    if angle_a < angle_b
        line = createLine(hull(p_a, :), hull(ind_a2, :));
        width = distancePointLine(hull(p_b, :), line);
        p_a = mod(p_a, nV) + 1;
    
    else
        line = createLine(hull(p_b, :), hull(ind_b2, :));
        width = distancePointLine(hull(p_a, :), line);
        p_b = mod(p_b, nV) + 1;

    end
    
    % update points, maximum width and corresponding angle if needed
    if width > max_width
        p1  = hull(p_a,:);
        p2  = hull(p_b,:);      
%         max_width = width;
%         min_angle = rotated_angle;
    end

end
%% Sanity check for points neighbouring p1,p2
p1s = [hull(p_a-1,:); p1 ; hull(p_a+1,:)];
p2s = [hull(p_b+1,:); p2 ; hull(p_b-1,:)];
%% 
for i = 1:3
    for j = 1:3
    d(i,j) = pdist([p1s(i,:) ; p2s(j,:)],'euclidean');
    end
end
%%  
    maximum = max(max(d));
    [row,column]=find(d==maximum);
%     [tmp, ID] = max(d);
    p1  = p1s(row,:);
    p2  = p2s(column,:);
