numfiles = 20;
mydata = cell(1, numfiles);

for k = 1:numfiles
  myfilename = sprintf('profiles\profile%d.txt', k);
  mydata{k} = importdata(myfilename);
end
%% 
out_file_profile = sprintf('Bladegen/blade-profile.curve'); fileID1 = fopen(out_file_profile,'w');
out_file_hub     = sprintf('Bladegen/blade-hub.curve');     fileID2 = fopen(out_file_hub,'w');
out_file_shroud  = sprintf('Bladegen/blade-shroud.curve');  fileID3 = fopen(out_file_shroud,'w');
%% 

for i = 1:size(radius,2)
     
%% Create Surface #2: the cylinder

[z,y,x] = cylinder(radius(i),3.5*n_points);                                 %cylinder(radius, no. of sides)
DT = delaunayTriangulation(40.*x(:), y(:), z(:));
[Surface2.faces, Surface2.vertices] = freeBoundary(DT);
% hold on
% S=Surface2; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'b') 

%% Find intersection curve and plot it
[intersect12, Surf12] = SurfaceIntersection(Surface1, Surface2);            %find intersection curve
S=Surf12; 
%trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'EdgeColor', 'g', 'FaceColor', 'g')
axis equal
%% 

hold on  

profile_curve = delaunayTriangulation(Surf12.vertices(:,1), Surf12.vertices(:,2), Surf12.vertices(:,3));
%% Rotating and translating points

raw_curve = (profile_curve.Points)*rotx;
raw_curve(:,1) = 300 + raw_curve(:,1);

%% 

unwrapped = unwrap(raw_curve);
[uw_profile.LE uw_profile.TE] = get_LE_TE(unwrapped.points);
[uw_profile.curve uw_profile.LE uw_profile.TE uw_profile.S1 uw_profile.S2] = set_circular_sequence(unwrapped.points, uw_profile.LE, uw_profile.TE);
uw_profile.camber = get_camber(uw_profile.S1, uw_profile.S2, uw_profile.LE, uw_profile.TE);

profile.LE = wrap(uw_profile.LE, unwrapped.r);
profile.TE = wrap(uw_profile.TE, unwrapped.r);
profile.S1 = wrap(uw_profile.S1, unwrapped.r);
profile.S2 = wrap(uw_profile.S2, unwrapped.r);
profile.camber = wrap(uw_profile.camber, unwrapped.r);
profile.curve = wrap(uw_profile.curve, unwrapped.r);
%% Plotting points and curves

% plot3(profile.curve(:,1),profile.curve(:,2),profile.curve(:,3),'.g');
% plot3(profile.camber(:,1),profile.camber(:,2),profile.camber(:,3),'.r');
% plot3(profile.LE(1), profile.LE(2), profile.LE(3), '*b');
% plot3(profile.TE(1), profile.TE(2), profile.TE(3), '*b');
% text(profile.LE(1), profile.LE(2), profile.LE(3),'LE','FontSize',8)
% text(profile.TE(1), profile.TE(2), profile.TE(3),'TE','FontSize',8)

%% Find % span and output coordinates in BladeGen format

%fig = figure;
span = 100*(radius(i)-input.r_hub)/(input.r_tip-input.r_hub)
profile.curve(:,1) = 300 + profile.curve(:,1);
p_out = (profile.curve)*rotx;
plot3(p_out(:,1),p_out(:,2),p_out(:,3),'.g');
plot3(p_out(1,1),p_out(1,2),p_out(1,3),'*b');
fprintf(fileID1,'#Profile %d at %12.8f%%\n',i,span);
fprintf(fileID1,'%12.8f  %12.8f  %12.8f  \n',0.1*((profile.curve).'));


%% Hub and tip camber lines output

if      i==1
        hub = profile;
elseif  i==size(radius,2)
        tip = profile;
end 

%% Creating blade-hub and blade-shroud curve files
% 
if      i==1
        
        blade_hub_start = linspace(x_domain_start,hub.LE(1),10) ;
        blade_hub_start(2,:) = hub.LE(2);
        blade_hub_start(3,:) = hub.LE(3);
        
        blade_hub_end = linspace(hub.TE(1),x_domain_end,10);
        blade_hub_end(2,:) = hub.TE(2);
        blade_hub_end(3,:) = hub.TE(3);
        
        blade_hub = [ blade_hub_start hub.camber.' blade_hub_end ]; 
        
        blade_hub(1,:) = 300 + blade_hub(1,:);
        fprintf(fileID2,' %12.8f  %12.8f  %12.8f  \n', 0.1*blade_hub );
        
elseif  i==size(radius,2)
        
        A = linspace(x_domain_start,tip.LE(1),10);
        B = linspace(tip.TE(1),x_domain_end,10);
        
        blade_tip_start = A ;
        blade_tip_start(2,:) = tip.LE(2);
        blade_tip_start(3,:) = tip.LE(3);
        
        blade_tip_end = B ;
        blade_tip_end(2,:) = tip.TE(2);
        blade_tip_end(3,:) = tip.TE(3);
        
        blade_tip = [blade_tip_start tip.camber.' blade_tip_end ]; 
        
        blade_tip(1,:) = 300 + blade_tip(1,:);
        fprintf(fileID3,' %12.8f  %12.8f  %12.8f  \n', 0.1*blade_tip );
end
%% 
end
%% 

% xlim([-100 200])
% ylim([-100 400])
% zlim([-150 150])

fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
%% 
% scatter3(blade_hub(1,:),blade_hub(2,:),blade_hub(3,:))
% scatter3(blade_tip(1,:),blade_tip(2,:),blade_tip(3,:))