function unwrapped = unwrap(raw_curve)
%% This function unwraps a constant radii 3D curve to a planar 2D curve.
%  
%  Input - 3D raw_curve with columns x,y,z
%  Output - points: 2D unwrapped curve with columns x, r*theta
%                r: radius of the curve          
%% 
vertices = raw_curve;
x = vertices(:,1);
y = vertices(:,2);
z = vertices(:,3);
%% 
[theta,rho,x] = cart2pol(y,z,x);
%% 
unwrapped.r = mean(rho);

unwrapped.points = [ x (unwrapped.r).*theta ];

%scatter(unwrapped.points(:,1), unwrapped.points(:,2))

end