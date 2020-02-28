function [camber_points, camber_normals] = create_camber(US,LS)
%% Calculate camber line and camber normals of 2D Whittle fan
% 1. Calculate halfway point between upper and lower surface for each
% x-location.
% 2. Add to an array with camber mean line points, with LE and TE
% points first and last respectively.
% 3. Interpolate line through all mean camber line points
% 4. Create camber normal values at every point using mean camber line
% slope at each x-location
%% 
% Step 1 and 2
LE = US(1,:);
TE = US(end,:);
camber_points = zeros(length(US),2);
camber_points(1,:) = LE;
camber_points(end,:) = TE;

for i = 2:(length(US)-1)
    x = US(i,1);
    y_us = US(i,2);
    y_ls = LS(i,2);
    delta = y_us-y_ls;
    y = y_ls + (delta/2);
    camber_points(i,1) = x;
    camber_points(i,2) = y;
end

% % Step 3
% num_points = 1000;
% x_start = US(1,1);
% x_end = US(end,1);
% interp_locations = linspace(x_start, x_end, num_points);
% camber_points_interp = interp1(camber_points(:,1), camber_points(:,2), interp_locations, 'spline');

% Step 4
camber_normals = zeros(length(US),3);
for i = 1:length(US)
   if i == 1
       dy = camber_points(i+1,2) - camber_points(i,2);
       dx = camber_points(i+1,1) - camber_points(i,1);
       angle = atan(dy/dx);
       Nx = -sin(angle);
       Ny = cos(angle);
       camber_normals(i,1) = camber_points(i,1);
       camber_normals(i,2) = Nx;
       camber_normals(i,3) = Ny;
   elseif i == length(US)
       dy = camber_points(i,2) - camber_points(i-1,2);
       dx = camber_points(i,1) - camber_points(i-1,1);
       angle = atan(dy/dx);
       Nx = -sin(angle);
       Ny = cos(angle);
       camber_normals(i,1) = camber_points(i,1);
       camber_normals(i,2) = Nx;
       camber_normals(i,3) = Ny;
   else
       dy = camber_points(i+1,2) - camber_points(i-1,2);
       dx = camber_points(i+1,1) - camber_points(i-1,1);
       angle = atan(dy/dx);
       Nx = -sin(angle);
       Ny = cos(angle);
       camber_normals(i,1) = camber_points(i,1);
       camber_normals(i,2) = Nx;
       camber_normals(i,3) = Ny;
   end
end

end

