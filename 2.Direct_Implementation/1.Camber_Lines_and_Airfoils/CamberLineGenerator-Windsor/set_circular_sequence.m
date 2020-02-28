function [LE TE S1 S2] = set_circular_sequence(xy_points, LE, TE)
% This function takes a curve (delaunay triangulation format), its approximate 
% leading and trailing edge as the input.
% 
% The function outputs the exact LE, TE, and the suction and pressure
% sides as different surfaces. It also outputs the curve as a set of points with
% a circular sequence.

% clearvars -except unwrapped LE TE
% xy_points = unwrapped.points;
%% 

x = xy_points(:,1);
y = xy_points(:,2);
%plot3(LE(:,1),LE(:,2), LE(:,3),'*r')
% hold on
%plot3(TE(:,1),TE(:,2), TE(:,3),'*r')
%% 
cx = mean(x);
cy = mean(y);
a = -atan2(y - cy, x - cx);
[~, order] = sort(a);
X = x(order);
Y = y(order);
%% 
airfoil= [X Y];
[p,ID1] = ismember(LE,airfoil,'rows');
[p,ID2] = ismember(TE,airfoil,'rows');
%% 
if ID1 < ID2
    
elseif ID2 < ID1
    
   temp = ID2;
   ID2 = ID1;
   ID1 = temp;
       
end
%% 

 S1 = airfoil(ID1+1:ID2-1,:);

%plot3(S1(:,1),S1(:,2),S1(:,3),'*b');
%axis equal
S2 = [ airfoil(ID2+1:end, :) ; airfoil(1:ID1-1, :) ];
%hold on
%plot3(S2(:,1),S2(:,2),S2(:,3),'*g');

new_curve = [ LE ; S1 ; TE ; S2 ; LE ];
end