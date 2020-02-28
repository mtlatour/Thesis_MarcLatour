clear 
clc

filename = fullfile('blade-hub.txt');
% filename = fullfile('blade-shroud.txt');
data = readtable(filename,'Format','%f%f%f','ReadVariableNames',false);
T = table2array(data);
%% 

x = T(:,1);
y = T(:,2);
z = T(:,3);
r = sqrt(x.^2 + y.^2 + z.^2);

for i=1:size(T,1);
    hold on
    plot3(T(i,1),T(i,2),T(i,3),'.b')
    scatter3(T(i,1),T(i,2),T(i,3),'og')
    pause(0.05);
end
