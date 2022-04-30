close all; clear all; clc; addpath('src');
%% Initialization
IX1 = 3; 
IY1 = 4; 
IZ1 = 38;
IT = 12;
IP = 43;

err = 1e-3;

%% Target

load('data/vind.mat');
load('data/x_range.mat');
load('data/y_range.mat');
load('data/z_range.mat');
load('data/theta_range.mat');
load('data/phi_range.mat');

X1 = x_range( IX1 );    Y1 = y_range( IY1 );    Z1 = z_range( IZ1 );

v1 = vind(IX1, IY1, IZ1, IT, IP);

%% Prediction
idx_preds = find( abs(vind-v1) < err);

%% Plot
figure

[ix1, iy1, iz1, ~, ~] = ind2sub(size(vind), idx_preds);
scatter3(x_range(ix1), y_range(iy1), z_range(iz1), '.');

hold on,
scatter3(X1, Y1, Z1, 'ro');
