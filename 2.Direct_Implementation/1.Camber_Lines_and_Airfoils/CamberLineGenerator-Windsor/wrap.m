function wrapped_points = wrap(points, radius)
%% This function unwraps a planar 2D curve to a constant radius 3D curve.
%  
%  Input - 2D unwrapped curve with columns x, r*theta
%          'radius' of the curve
%  Output - 3D raw_curve with columns x,y,z 
%% 
x = points(:,1);
rtheta = points(:,2);
rho    = repmat(radius,[size(x),1]);
theta  = rtheta./radius;

%% 
[y,z,x] = pol2cart(theta,rho,x);
%% 
wrapped_points = [ x y z ];

% plot3(x,y,z,'.r')
%scatter(unwrapped_points(:,1), unwrapped_points(:,2))

end