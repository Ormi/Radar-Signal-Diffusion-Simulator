% // M = [0 0 0.5; 
% //      0 -1 0; 
% //      2 -1 0;
% //      0 1 0.5; 
% //      0 2 0; 
% //      0 -2 0; 
% //      2 0 0.5; 
% //      2 -2 0; 
% //      2 1 0.5; 
% //      2 2 0;
% //      2 -2 0; 
% //      1 -2 0; 
% //      1 2 0;
% //      1 0 0.5;
% //      1 1 0.5];

clear all;
addpath(genpath('./jsonlab'))
savepath

%%
% Setting static variables for simulation
data_model=loadjson('model_bike.json');

M = (data_model.directions);
x = M(:,1);
y = M(:,2);
z = M(:,3);

xlin = linspace(min(x),max(x),33);
ylin = linspace(min(y),max(y),33);

[X,Y] = meshgrid(xlin,ylin);

f = scatteredInterpolant(x,y,z);
Z = f(X,Y);

figure
mesh(X,Y,Z) %interpolated
axis tight; hold on
plot3(x,y,z,'.','MarkerSize',35) %nonuniform