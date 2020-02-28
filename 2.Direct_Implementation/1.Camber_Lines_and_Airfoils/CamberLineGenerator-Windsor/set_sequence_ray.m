function [LE TE S1 S2] = set_sequence_ray(xy_points, LE, TE)
%% 
% This function puts the profile points in a circular sequence using a
% revolving ray originating from a point defined by the mean x-y coordinates.
% 
% Inputs -  2D profile curve with two columns x,y 
%           Leading edge of airfoil profile
%           Trailing edge of airfoil profile
% 
% Outputs - Same LE, TE
%           Suction/Pressure sides S1,S2 as seperate sequenced curves with x,y columns
% 
% This function might give faulty results with complex profile shapes.
% Use function correct_sequence to reorder points based on the closest 
% neighbour (minimum distance) method.
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

% new_curve = [ LE ; S1 ; TE ; S2 ; LE ];
end