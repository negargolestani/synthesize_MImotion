close all; clear all; clc; addpath('src');
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
folder_path = strcat('../datasets/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));
file_name = file_list( 1 ).name;

data = readtable( strcat(folder_path, file_name) );
% data = data(50:end,:);
Nt =  height(data);

for t = 1:Nt
    C1(t,:) = str2num(erase(data.center_1{t}, ["[","]"]))*100;
    norm(t,:) = str2num(erase(data.norm{t}, ["[","]"]));
    theta(t) = acos(norm(t, 3));
    phi(t) = mod(atan2(norm(t,2),norm(t,1)), 2*pi);
end

load('data/vind.mat');
load('data/x_range.mat');
load('data/y_range.mat');
load('data/z_range.mat');
load('data/theta_range.mat');
load('data/phi_range.mat');
[X, Y, Z, THETA, PHI] = ndgrid(x_range, y_range, z_range, theta_range, phi_range);
DIST = sqrt( X.^2 + Y.^2 + Z.^2 );

err = .001;
err_dist = 1;

%% Target
ix_trg=18; iy_trg=4; iz_trg=14; it_trg=5; ip_trg=26;

x_trg = x_range(ix_trg); 
y_trg = y_range(iy_trg); 
z_trg = z_range(iz_trg);

theta_trg = theta_range(it_trg);  
phi_trg = phi_range(ip_trg);    
n_trg = [sin(theta_trg)*cos(phi_trg), sin(theta_trg)*sin(phi_trg), cos(theta_trg)];

c1_trg = [x_trg, y_trg, z_trg];
c2_trg = c1_trg + 3 * n_trg;
    
V1 = round( vind(ix_trg, iy_trg, iz_trg, it_trg, ip_trg), 3);
V2 = round( synthesize_vind(c2_trg*1e-2, n_trg), 3);

d_est = sqrt( x_trg^2 + y_trg^2 + z_trg^2 );
%%
cond1 = abs(DIST-d_est) < err_dist;
cond2 = vind~=0 & abs(vind-V1)< err;
idx_list = find( cond1 .* cond2 );

sol = [];
for i = 1:length(idx_list)
    [ix, iy, iz, it, ip] = ind2sub(size(vind), idx_list(i));
    x = x_range(ix); y = y_range(iy); z = z_range(iz);
    theta = theta_range(it);  phi = phi_range(ip);
    
    n = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    c1 = [x, y, z] ;
    c2 = c1 + 3 * n;
    v2 = synthesize_vind(c2*1e-2, n);
    
    if abs(V2 - v2)< err
        sol = [sol; x, y, z, theta, phi];
    end
end

%%
figure, 
scatter3(sol(:,1), sol(:,2), sol(:,3), 'b.')
hold on, 
scatter3(c1_trg(1), c1_trg(2), c1_trg(3), 'ro')


