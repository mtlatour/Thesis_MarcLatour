tic
clear
% clearvars -except allprofiles
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_points = 500;                                                             % no. of points for each airfoil curve
dt = stlread('Pacifica_single_blade_mesh.stl');                             % .stl file to be read

%% dimensions are in millimetres
input.r_hub = 88;                                                         % approximate hub radius
input.r1 = 92;                                                              % filleting end at hub
input.r2 = 222;                                                             % filleting start at tip
input.r_tip = 227;                                                          % approximate tip radius
n_profiles.hub = 7;                                                         % no. of profiles b/w r_hub and r1
n_profiles.mid = 21;                                                        % no. of profiles b/w r1 and r2
n_profiles.tip = 8;                                                         % no. of profiles b/w r2 and r_tip
%% Create Surface #1: the fan
Surface1.faces = [dt.ConnectivityList]; Surface1.vertices = [dt.Points(:,1),dt.Points(:,2),dt.Points(:,3)];
% S=Surface1; trisurf(S.faces,0.1*S.vertices(:,1),0.1*S.vertices(:,3),0.1*S.vertices(:,2),'FaceAlpha', 0.2, 'FaceColor', 'b')
axis equal
%% Create radii vector
nearhub_radii = linspace(input.r_hub, input.r1, n_profiles.hub);
mid_radii     = linspace(input.r1, input.r2, n_profiles.mid); 
neartip_radii = linspace(input.r2, input.r_tip, n_profiles.tip);
radius        = [nearhub_radii mid_radii neartip_radii];                    %Concatenating vectors
radius = unique(radius);
%% 
for i = 1:size(radius,2)
%% Create Surface #2: the cylinder
[z,y,x] = cylinder(radius(i),3.5*n_points);                                 %cylinder(radius, no. of sides)
DT = delaunayTriangulation(40.*x(:), y(:), z(:));
[Surface2.faces, Surface2.vertices] = freeBoundary(DT);
% S=Surface2; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'b') 
%% Find intersection curve and plot it
[intersect12, Surf12] = SurfaceIntersection(Surface1, Surface2);            %find intersection curve
S=Surf12; %trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'EdgeColor', 'g', 'FaceColor', 'g')
profile_curve = delaunayTriangulation(Surf12.vertices(:,1), Surf12.vertices(:,2), Surf12.vertices(:,3));
raw_curve = (profile_curve.Points);
%% Unwrap and process curve

%converts (x,y,z) to (x,r*theta)
unwrapped = unwrap(raw_curve);       
%find LE,TE
[uw_profile.LE uw_profile.TE] = max_distance_callipers(unwrapped.points);  
%reorder profile points in circular sequence using a revolving ray
[uw_profile.LE uw_profile.TE uw_profile.S1 uw_profile.S2] = set_sequence_ray(unwrapped.points, uw_profile.LE, uw_profile.TE);
%correct regions with faulty sequence using the nearest neighbour method
[uw_profile.S1] = correct_sequence(uw_profile.S1, uw_profile.LE);
uw_profile.curve = [uw_profile.LE ; uw_profile.S1 ; uw_profile.TE ; uw_profile.S2 ; uw_profile.LE]; 
%smooth uneven leading and trailing edge region
uw_profile.curve = smooth_curve (uw_profile.curve, uw_profile.TE, 'TE', 30);
uw_profile.curve = smooth_curve (uw_profile.curve, uw_profile.LE, 'LE', 80);
uw_profile.S2 = smooth_curve (uw_profile.S2, uw_profile.LE, 'LE', 80);
%get camber line
uw_profile.camber = get_accurate_camber(uw_profile.S1, uw_profile.S2, uw_profile.LE, uw_profile.TE, 'linear', 50);
% %%  
% fig = figure
% hold on; axis equal; axis tight;  grid on
% plot(uw_profile.curve(:,1).',uw_profile.curve(:,2).','g');
% % rainbowplot(uw_profile.curve(:,1).',uw_profile.curve(:,2).');
% plot(uw_profile.camber(:,1),uw_profile.camber(:,2),'.r');
% plot(uw_profile.LE(1), uw_profile.LE(2),  'ok',uw_profile.TE(1), uw_profile.TE(2), 'ok');

%% Wrap back all points and convert everything to centimetres.
profile.LE = 0.1*wrap(uw_profile.LE, unwrapped.r);          profile.TE = 0.1*wrap(uw_profile.TE, unwrapped.r);
profile.S1 = 0.1*wrap(uw_profile.S1, unwrapped.r);          profile.S2 = 0.1*wrap(uw_profile.S2, unwrapped.r);
profile.camber = 0.1*wrap(uw_profile.camber, unwrapped.r);  profile.curve = 0.1*wrap(uw_profile.curve, unwrapped.r);
%% 
% hold on; axis equal; grid on
% fig = figure
% plot3(profile.curve(:,1).',profile.curve(:,3).',profile.curve(:,2).','--b','LineWidth',2);
% rainbowplot(uw_profile.curve(:,1).',uw_profile.curve(:,2).');
% plot(uw_profile.camber(:,1),uw_profile.camber(:,2),'.r');
% plot(uw_profile.LE(1), uw_profile.LE(2),  'ok',uw_profile.TE(1), uw_profile.TE(2), 'ok');
%% Find percent span of profile.
span = 100*(radius(i)-input.r_hub)/(input.r_tip-input.r_hub)
%% Create a cell array with profile information - i)radius ii)span iii)profile details
allprofiles.radius{i}   = unwrapped.r;
allprofiles.span{i}     = span;
allprofiles.profile{i}  = profile;
uwprofiles{i}.profile{i} = uw_profile;
end
%% 
toc
save ('allprofiles.mat', 'allprofiles')

save ('uwprofiles.mat', 'uwprofiles')

