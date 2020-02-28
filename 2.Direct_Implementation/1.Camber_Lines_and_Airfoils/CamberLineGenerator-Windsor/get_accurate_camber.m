function camber = get_accurate_camber(US, LS, point_LE, point_TE, type, n)
%% This function estimates the camber line of an airfoil profile.
% It uses an iterative method to estimate a camber line having equal 
% perpendicular distances from the pressure and suction surface at all
% points. The points can have a linear or cosine spacing (denser near LE/TE).
% 
% Input - Suction/Pressure 2D surfaces US, LS with x,y columns
%         Leading and Trailing edge points LE, TE
%         'type' can be 'linear' or 'cosine'
%         No. of points 'n'
%         
% Output - Camber line with x,y columns
% 
% Example : 
% camber = get_accurate_camber(pressure_surf,suction_surf,LE,TE,'linear',60)
% 
% This function is dependent on the function get_approx_camber.
% Use function max_distance_callipers to find the LE-TE.
%% 
approx_cam = get_approx_camber(US, LS, point_LE, point_TE, type, n);
LE = point_LE'; 
TE = point_TE';
S1 = [LE  US'  TE];
S2 = [TE  US'  LE];
%% 
n = 2:(size(approx_cam,2)-1);
n = flip(n);
iter=5;
tol = 10;
cam = approx_cam';
cam_new = approx_cam';
%% 
while (tol>=0.001)

for i = 1:size(n,2)
    k=n(i);
if k~=0
count1 = 0;
count2 = 0;
p1 = cam_new(:,k-1);
p2 = cam_new(:,k+1);
pm = cam_new(:,k);

%% Drawing line perpendicular to p1-p2
m = - (p1(1)-p2(1))/(p1(2) - p2(2));
b = pm(2) - m*pm(1);
%% Find interesection with curves
y_points = linspace(m*(pm(1)-5)+b,m*(pm(1)+5)+b,10);
x_points = (y_points - b)/m;
L = [x_points ; y_points];
%% 

P1=InterX(L,S1);
P2=InterX(L,S2);
PM=0.5*(P1+P2);

if          size(PM,2)==0;
                PC = cam_new(:,k);
elseif      size(PM,2)==1     
               	PC = PM;
elseif      size(PM,2)>1
                diff = PM - repmat(pm, 1, length(PM));
                d = diff(1,:).^2 + diff (2,:).^2;
                [M,I] = min(d);
                PC = PM(:,I);
end

cam_new(:,k) = PC;
end
end

    iter = iter + 1;
    
    cam_diff = cam_new - cam;
    tol = max(cam_diff(1,:).^2 + cam_diff (2,:).^2);
    cam = cam_new;
end

camber = cam';

end