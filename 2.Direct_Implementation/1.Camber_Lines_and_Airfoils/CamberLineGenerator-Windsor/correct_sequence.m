function new_S = correct_sequence (S , LE)
%% 
% This function puts the profile points in an ordered sequence based on the 
% closest neighbour (minimum distance) method.
% 
% Inputs -  2D profile curve with two columns x,y 
%           Desired first point on the curve
% 
% Outputs - Reordered 2D profile curve with two columns x,y 
% 
% This method, though more accurate, is slower than the set_sequence_ray
% function which performs the same job but can give faulty results 
% with complex profile shapes.
%% 
x = S(:,1).';
y = S(:,2).';

% rainbowplot(x,y)

temp_x = LE(1);
temp_y = LE(2);

X = temp_x;
Y = temp_y;
%% 
for i=1:size(x,2)
%% 
for j = 1:size(x,2) 
        x_diff(j) = temp_x - x(j);            
        y_diff(j) = temp_y - y(j); 
        d = x_diff.^2 + y_diff.^2;
    [M1, ID] = min(d(d>0));
end
%% 
                                                        %max value in each row
    temp_x = x(ID);
    temp_y = y(ID);
     X = [X temp_x];
     Y = [Y temp_y];
     x(ID) = [];
     y(ID) = [];
end
%% 

% rainbowplot(X,Y)

new_S = [X(2:end) ; Y(2:end)].';
end