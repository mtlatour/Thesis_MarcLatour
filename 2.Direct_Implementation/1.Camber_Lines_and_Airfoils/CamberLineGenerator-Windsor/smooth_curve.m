function curve = smooth_curve (curve , point, type, n)
%% 
% This function creates a convex hull around the neighbouring region of the 
% specified point on a curve by fitting a convex hull on it.
% 
% Inputs -  'curve' : 2D profile curve with two columns x,y 
%           'point' : coordinates of the specified point
%           'type'  : 'LE' or 'TE'
%           'n' : no. of points on each side of the specified point that
%           are to be included in computing the convex hull.

% Outputs - 'curve' : Reconstructed 2D profile curve with two columns x,y 
%% 
% clear 
% clc
% load uw_profile.mat
% 
% curve = uw_profile.curve;
% point = uw_profile.TE;
% n=10;
% type = 'LE';
%% find points of interest
if (type=='TE')
I = find(curve==point);
c_points = curve(I(1)-n:I(1)+n, :);
elseif (type=='LE')
c_points = [curve(end-n:end-1, :) ; curve(1:n, :)];
end
%% find convex hull of points of interest
inds = convhull(c_points(:,1),c_points(:,2));                               %indices of points in 'c_points' to keep
hull = c_points(inds, :);                                                   %points of interest to keep
%% find points which do not lie on the convex hull
temp = find(~ismember(c_points,hull,'rows'));
bad_points = c_points(temp, :);
%% delete bad points from the curve
bad_indices = find(ismember(curve,bad_points,'rows'));
curve(bad_indices, :) = [];
end