clear
clc
load allprofiles.mat

x_domain_start = -300;                                                      %inlet  x-plane
x_domain_end   = 500;                                                       %outlet x-plane

%% Writing all profile curve points sequenced from lE to TE to LE in Bladegen format

out_file_profile = sprintf('Bladegen/blade-profile.curve'); fileID = fopen(out_file_profile,'w');
for i = 1:size(allprofiles.profile,2)
fprintf(fileID,'#Profile %d at %12.8f%%\n' , i , allprofiles.span{i});
fprintf(fileID,'%12.8f  %12.8f  %12.8f \n',  allprofiles.profile{i}.curve.');
end
fclose(fileID);

%% Writing hub-curve file for Bladegen.

hub = allprofiles.profile{1};                                               % hub profile details
hub_upstream = linspace(0.1*x_domain_start,hub.LE(1),10) ;                      % x-coord of upstream hub line
hub_upstream(2,:) = hub.LE(2); hub_upstream(3,:) = hub.LE(3);               % y,z-coord of downstream hub line
hub_cam = hub.camber.';                                     % hub blade camber line
hub_downstream = linspace(hub.TE(1),0.1*x_domain_end,10);                       % x-coord of downstream hub line
hub_downstream(2,:) = hub.TE(2);   hub_downstream(3,:) = hub.TE(3);           % y,z-coord of downstream hub line

hub_line = [ hub_upstream hub_cam hub_downstream ]; 

out_file_hub = sprintf('Bladegen/blade-hub.curve');     fileID = fopen(out_file_hub,'w');
fprintf(fileID,' %12.8f  %12.8f  %12.8f  \n', hub_line );
fclose(fileID);
%% Writing shroud-curve file for Bladegen.

tip = allprofiles.profile{end};                                               % shroud profile details
tip_upstream = linspace(0.1*x_domain_start,tip.LE(1),10) ;                      % x-coord of upstream shroud line
tip_upstream(2,:) = tip.LE(2); tip_upstream(3,:) = tip.LE(3);               % y,z-coord of downstream shroud line
tip_cam = tip.camber.';                                     % shroud blade camber line
tip_downstream = linspace(tip.TE(1),0.1*x_domain_end,10);                       % x-coord of downstream shroud line
tip_downstream(2,:) = tip.TE(2);   tip_downstream(3,:) = tip.TE(3);           % y,z-coord of downstream shroud line

tip_line = [ tip_upstream tip_cam tip_downstream ]; 

out_file_shroud  = sprintf('Bladegen/blade-shroud.curve');  fileID = fopen(out_file_shroud,'w');
fprintf(fileID,' %12.8f  %12.8f  %12.8f  \n', tip_line );
fclose(fileID);
%% Writing camber files
% 
% fprintf(fileID,' %12.8f  %12.8f  %12.8f  \n', uw_profile.camber' );
% fclose(fileID);
% 
% out_file = sprintf('camber/profile%d',i); fileID = fopen(out_file,'w');
%% Plotting 

for i=1:size(allprofiles.profile,2)
    
    S = allprofiles.profile{i}.curve.';
    plot3(S(1,:),S(2,:),S(3,:), 'g');
    hold on
    axis equal
    view(3)
    grid on
end
%% 

  plot3(hub_line(1,:),hub_line(2,:),hub_line(3,:), 'k')
  plot3(tip_line(1,:),tip_line(2,:),tip_line(3,:), 'b');
%% 
    hold on
    axis equal
    view(3)
    grid on
plot3(hub_cam(1,:),hub_cam(2,:),hub_cam(3,:), 'r');
plot3(hub_upstream(1,:),hub_upstream(2,:),hub_upstream(3,:), 'r');
plot3(hub_downstream(1,:),hub_downstream(2,:),hub_downstream(3,:), 'r');
