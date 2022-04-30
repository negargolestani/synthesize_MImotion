close all; clear all; clc; addpath('src');
%% Initialization
IX1 = 3; 
IY1 = 4; 
IZ1 = 38;
IT = 12;
IP = 43;

err = 1e-3;
D=3;

%% Target

load('data/vind.mat');
load('data/x_range.mat');
load('data/y_range.mat');
load('data/z_range.mat');
load('data/theta_range.mat');
load('data/phi_range.mat');

X1 = x_range( IX1 );    Y1 = y_range( IY1 );    Z1 = z_range( IZ1 );
THETA = theta_range( IT );   PHI = phi_range( IP );
C2 = [X1, Y1, Z1] + 3*[sin(THETA)*cos(PHI), sin(THETA)*sin(PHI), cos(THETA)];
[~,IX2] = min(abs(C2(1)-x_range)); 
[~,IY2] = min(abs(C2(2)-y_range)); 
[~,IZ2] = min(abs(C2(3)-z_range)); 

v1 = vind(IX1, IY1, IZ1, IT, IP);
v2 = vind(IX2, IY2, IZ2, IT, IP);       

%% Prediction
idx_preds_1 = find( abs(vind-v1) < err);
idx_preds = [];

for i = 1:length(idx_preds_1)
    
    [ix1, iy1, iz1, it, ip] = ind2sub(size(vind), idx_preds_1(i));
    theta = theta_range(it);   phi = phi_range(ip);
    c2 = [x_range(ix1), y_range(iy1), z_range(iz1)] + ...
        D*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];     
    [~,ix2] = min(abs(c2(1)-x_range)); 
    [~,iy2] = min(abs(c2(2)-y_range)); 
    [~,iz2] = min(abs(c2(3)-z_range)); 
    
    if abs(vind(ix2, iy2, iz2, it, ip)-v2) < err
        idx_preds = [idx_preds; idx_preds_1(i)];
    end
    
end
  
%% Plot
figure

[ix1, iy1, iz1, ~, ~] = ind2sub(size(vind), idx_preds);
scatter3(x_range(ix1), y_range(iy1), z_range(iz1), '.');

hold on,
scatter3(X1, Y1, Z1, 'ro');
