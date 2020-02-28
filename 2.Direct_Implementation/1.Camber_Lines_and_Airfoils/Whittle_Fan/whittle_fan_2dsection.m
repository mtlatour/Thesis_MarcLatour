% Code to create a 2D section of the whittle fan at mid-span
clear;
clc;
load('bladeshapes.mat')
set(0,'defaulttextinterpreter','latex')

% hold on 
% for i = 1:54
%     PSx = RPSx(:,i); PSy = RPSy(:,i); PSz = RPSz(:,i); SSx = RSSx(:,i);
%     SSy = RSSy(:,i); SSz = RSSz(:,i); if i == 27
%         scatter3(PSx, PSy, PSz, 'b*') scatter3(SSx, SSy, SSz, 'b*')
%     else
%         scatter3(PSx, PSy, PSz, 'k*') scatter3(SSx, SSy, SSz, 'r*')
%     end
% end

% Middle span-fraction 3D section
midspan = 27;
mid_ps_x = RPSx(:,midspan);
mid_ps_y = RPSy(:,midspan);
mid_ps_z = RPSz(:,midspan);
mid_ss_x = RSSx(:,midspan);
mid_ss_y = RSSy(:,midspan);
mid_ss_z = RSSz(:,midspan);

% % 3D scatter of mid-span section
% hold on
% scatter3(mid_ps_x, mid_ps_y, mid_ps_z, 'k*')
% scatter3(mid_ss_x, mid_ss_y, mid_ss_z, 'r*')
% hold off

% 2D plot of mid-span section
min_y_ps = max(mid_ps_y);
min_y_ss = max(mid_ss_y);
min_y = max(min_y_ps,min_y_ss);
min_x = min(mid_ps_x);
figure()
hold on
plot(mid_ps_x - min_x, -mid_ps_y + min_y, 'k*')
plot(mid_ss_x - min_x, -mid_ss_y + min_y, 'r*')
xlabel("x")
ylabel("y")

% Create correctly oriented x,y coordinates of both surfaces and save to
% txt files
uppersurface = [mid_ss_x - min_x, -mid_ss_y + min_y];
lowersurface = [mid_ps_x - min_x, -mid_ps_y + min_y];
uppersurface_alignedpoints = uppersurface;
uppersurface_alignedpoints(1,:) = lowersurface(1,:);
uppersurface_alignedpoints(end,:) = lowersurface(end,:);
save('Whittle_upper.txt', 'uppersurface', '-ascii')
save('Whittle_lower.txt', 'lowersurface', '-ascii')
save('Whittle_upper_alignedendpoints.txt', 'uppersurface_alignedpoints', '-ascii')

% Get camber points and normals of 2D section
[camberpoints, cambernormals] = create_camber(uppersurface_alignedpoints, lowersurface);
figure()
hold on
upper = plot(uppersurface_alignedpoints(:,1),uppersurface_alignedpoints(:,2), 'k-', 'linewidth', 2.5);
lower = plot(lowersurface(:,1), lowersurface(:,2), 'k-', 'linewidth', 2.5);
camber = plot(camberpoints(:,1), camberpoints(:,2), 'r-', 'linewidth', 2.5);
xlabel("x")
ylabel("y")
% title("2D Section of Whittle Fan, taken at constant span fraction")
grid()
set(gca, 'fontsize', 14)
legend([upper, camber], {'Airfoil','Camber'}, 'location', 'northwest')

save('Whittle_2D_cambernormals.txt', 'cambernormals', '-ascii')

% Shift coordinates for use in physical blade simulation
pitch = 0.05105;
upper_y_shifted = uppersurface_alignedpoints(:,2) + pitch/2;
lower_y_shifted = lowersurface(:,2) + pitch/2;
upper_shifted = [uppersurface_alignedpoints(:,1), upper_y_shifted];
lower_shifted = [lowersurface(:,1), lower_y_shifted];
save('whittle_upper_shifted.txt', 'upper_shifted', '-ascii')
save('whittle_lower_shifted.txt', 'lower_shifted', '-ascii')
