close all; clear all; clc; addpath('src');
%% Initialization
IX1 = 3;
IY1 = 4;
IZ1 = 38;
IT1 = 12;
IP1 = 43;
IP2 = 12;

err = 1e-3;
D = 3;

%% Target

load('data/vind.mat');
load('data/x_range.mat');
load('data/y_range.mat');
load('data/z_range.mat');
load('data/theta_range.mat');
load('data/phi_range.mat');
[X, Y, Z, THETA, PHI] = ndgrid(x_range, y_range, z_range, theta_range, phi_range);

X1 = x_range( IX1 );    Y1 = y_range( IY1 );    Z1 = z_range( IZ1 );
THETA1 = theta_range( IT1 );   PHI1 = phi_range( IP1 );
PHI2 = phi_range(IP2);  % can be any value
[~,IT2] = min(abs(acot(-cos(PHI1-PHI2)*tan(THETA1))-theta_range));
THETA2 = theta_range(IT2);
NORM2 = [sin(THETA2)*cos(PHI2), sin(THETA2)*sin(PHI2), cos(THETA2)];
C2 = [X1, Y1, Z1] - D * NORM2;
[~,IX2] = min(abs(C2(1)-x_range));
[~,IY2] = min(abs(C2(2)-y_range));
[~,IZ2] = min(abs(C2(3)-z_range));


v1 = vind(IX1, IY1, IZ1, IT1, IP1);
v2 = vind(IX2, IY2, IZ2, IT2, IP2);

%% Prediction
idx_preds_1 = find( abs(vind-v1) < err);
idx_preds = [];
for i = 1:length(idx_preds_1)
    
    [ix1, iy1, iz1, it1, ip1] = ind2sub(size(vind), idx_preds_1(i));
    x1 = x_range( ix1 );    y1 = y_range( iy1 );    z1 = z_range( iz1 );
    theta1 = theta_range(it1);   phi1 = phi_range(ip1);

    for ip2 = 1:length(phi_range)
        phi2 = phi_range(ip2);  % can be any value
        theta2 = acot(-cos(phi1-phi2)*tan(theta1));
        if theta2 < theta_range(1) || theta2 > theta_range(end)
            continue
        end
        [~,it2] = min(abs(theta2-theta_range));
        theta2 = theta_range(it2);
        norm2 = [sin(theta2)*cos(phi2), sin(theta2)*sin(phi2), cos(theta2)];
        c2 = [x1, y1, z1] - D * norm2;
        [~,ix2] = min(abs(c2(1)-x_range));
        [~,iy2] = min(abs(c2(2)-y_range));
        [~,iz2] = min(abs(c2(3)-z_range));
        
        
        if abs(vind(ix2, iy2, iz2, it2, ip2)-v2) < err
            idx_preds = [idx_preds; idx_preds_1(i)];
            break
        end
    end    
end

%% Plot
figure, set(gcf, 'Units', 'Inches', 'Position', [4,3,5,4]);

[ix1, iy1, iz1, ~, ~] = ind2sub(size(vind), idx_preds);
scatter3(x_range(ix1), y_range(iy1), z_range(iz1), '.');

hold on,
scatter3(X1, Y1, Z1, 'ro');
