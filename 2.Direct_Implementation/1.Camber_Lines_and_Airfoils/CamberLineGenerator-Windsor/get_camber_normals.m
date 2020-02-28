clear
clc
load allprofiles.mat

%% 
x = [];
y = [];
z = [];
LE = [];
TE = [];
%% 
figure
hold on


%% Constructing the camber line matrix
for i = 1:size(allprofiles.profile,2)
camber_line = allprofiles.profile{i}.camber;
x = [x , allprofiles.profile{i}.camber(2:end-1,1)];   % Constructing x-coord matrix
y = [y , allprofiles.profile{i}.camber(2:end-1,2)];   % Constructing y-coord matrix
z = [z , allprofiles.profile{i}.camber(2:end-1,3)];   % Constructing z-coord matrix 
LE = [LE ; allprofiles.profile{i}.LE];
TE = [TE ; allprofiles.profile{i}.TE];
hold on
end
% 
plot3(LE(:,1),LE(:,2),LE(:,3),'g','LineWidth',2)
hold on
plot3(TE(:,1),TE(:,2),TE(:,3),'r','LineWidth',2)

%% 
[ny, nz, nx] = surfnorm(y,z,x);
len = sqrt (nx.^2 + ny.^2 + nz.^2);
nx = nx./len;
ny = ny./len;
nz = nz./len; 
%% 
% plot surface and normals
surf(x,y,z)
hold on
quiver3(x,y,z,0.5*nx,0.5*ny,0.5*nz,'r')

%% Calculating nr and nth

nr  = ny.*y+nz.*z;
nth = nz.*y-ny.*z;
%% Reshaping matrices to a single row
X =reshape(x,1,[]);
Y =reshape(y,1,[]);
Z =reshape(z,1,[]);
R = sqrt(Y.^2+Z.^2);
NX =reshape(nx,1,[]);
NTH =reshape(nth,1,[]);
NR =reshape(nr,1,[]);

%% Construct openfoam-readable data files

files = ["nx","nth","nr"];
%% 
for i=1:3
    
    if i==1
        A = NX;
    elseif i==2
        A = NTH;
    elseif i==3
        A = NR;
    end
        
out_file = sprintf('%s_data', files(i));
fileID = fopen(out_file,'w');
fprintf(fileID,'//(x RADIUS %s ) \n ( \n', files(i));
fprintf(fileID,'( %12.8f  %12.8f  %12.8f ) \n', [X ; R ; A] );
fprintf(fileID, '%6s\n',')');
fclose(fileID);

end
