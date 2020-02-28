function [LE , TE] = get_LE_TE(xy_points)
% % % This function takes a curve (delaunay triangulation format), its approximate 
% % % leading and trailing edge as the input.
% % % The function outputs the exact LE, TE, and the suction and pressure
% % % sides as different surfaces. It also outputs the curve as a set of points with
% % % a circular sequence.
% curve     = profile_curve;
% xy_points = unwrapped_profile;
% clearvars -except unwrapped_profile;
%% 
x = xy_points(:,1);
y = xy_points(:,2);
%% 

for i=1:size(x,1)   
    for j=i:size(x,1)        
        x_diff(i,j) = x(i) - x(j);            
        y_diff(i,j) = y(i) - y(j); 
        d = x_diff.^2 + y_diff.^2;
    end   
end

[M1, ID] = max(d,[],2)                                                        %max value in each row

[M2,row] = max(M1);

column = ID(row);

%% 

x1 = x(row);
x2 = x(column);

if (x1<x2)
    LE = [x(row)    y(row)   ];
    TE = [x(column) y(column)];
else (x2<x1)
    LE = [x(column) y(column)];
    TE = [x(row)    y(row)   ];
end
    

% scatter(x,y,'g')
hold on
plot(LE(1),LE(2), 'or')
plot(TE(1),TE(2), '*b')
axis equal

end