clear
clc
load uwprofiles.mat

%% 

for i=1:size(uwprofiles,2)
x = uwprofiles{i}.profile{i}.camber(:,1);
y = uwprofiles{i}.profile{i}.camber(:,2);
m(i) = (y(2) - y(3))/(x(2) - x(3));
theta(i)= 90 - atand(m(i));
end
%% 

midtheta = theta(:,7:27);

