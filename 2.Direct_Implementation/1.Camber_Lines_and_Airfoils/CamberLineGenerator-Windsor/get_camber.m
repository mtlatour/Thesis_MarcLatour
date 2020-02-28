function camber = get_camber(US, LS, LE, TE)
%% 

S1 = [LE ; US ; TE];
S2 = [TE ; LS ; LE];

%% 

x1 = S1(:,1);
y1 = S1(:,2);

x2 = S2(:,1);
y2 = S2(:,2);

hold on
axis equal
% plot(x1.',y1.','.b');
% plot(x2.',y2.','.b');
%% 
v1 = [x1.'; y1.'];
v2 = [x2.'; y2.'];

% choose a point which will be the center of rotation
x_center = LE(1);
y_center = LE(2);


% choose angle of rotation
theta = atan2((LE(2)- TE(2)),(LE(1)- TE(1)));    

% create a matrix which will be used later in calculations
center1 = repmat([x_center; y_center], 1, length(x1));
center2 = repmat([x_center; y_center], 1, length(x2));

% define a counter-clockwise rotation matrix 
R1 = [cos(theta) sin(theta); -sin(theta) cos(theta)];

% do the rotation
s1  = v1 - center1;     % shift points in the plane so that the center of rotation is at the origin
s2  = v2 - center2;     % shift points in the plane so that the center of rotation is at the origin
so1 = R1*s1;             % apply the rotation about the origin
so2 = R1*s2;             % apply the rotation about the origin
% so1 = so1 + center1;
% so2 = so2 + center2;

x1_rotated = so1(1,:);
y1_rotated = so1(2,:);

x2_rotated = so2(1,:);
y2_rotated = so2(2,:);
%% 
% hold on
% 
% plot(LE(1),LE(2),'*g')
% plot(TE(1),TE(2),'*g')
% plot(x1_rotated, y1_rotated, 'r')
% plot(x2_rotated, y2_rotated, 'r')
%% Find camber surface

mn = min(min(x1_rotated(:,1)), min(x2_rotated(:,1)));
mx = max(max(x1_rotated(:,1)), max(x2_rotated(:,1)));
%query_points =  (mx - mn)*(0.5*(1-cos(linspace(0, pi, 32)))) + mn;
query_points =  linspace(mn,mx,60);

Y1_rotated = interp1(x1_rotated,y1_rotated,query_points);
Y2_rotated = interp1(x2_rotated,y2_rotated,query_points);

camber_rotated_X = query_points;
camber_rotated_Y = 0.5*(Y1_rotated + Y2_rotated);

% 
%% 

theta2 = -theta;

R2 = [cos(theta2) sin(theta2); -sin(theta2) cos(theta2)];

center3 = repmat([x_center; y_center], 1, length(camber_rotated_X));
v3 = [camber_rotated_X; camber_rotated_Y];
so3 = ((R2*v3)+ center3).'; 

camber_X = so3(:,1);
camber_Y = so3(:,2);

%% 
hold on
%plot(camber_rotated_X , camber_rotated_Y);

camber = [TE ; camber_X, camber_Y; LE];
camber = flip(camber);
scatter(camber(:,1) , camber(:,2), '.r');

end