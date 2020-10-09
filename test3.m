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

err = 2;
speed_max = 30;    % maximum velocity (cm/s)                               ref: https://watermark.silverchair.com/M584.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAr0wggK5BgkqhkiG9w0BBwagggKqMIICpgIBADCCAp8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMnRZ6g5KHto7PK_-kAgEQgIICcOeLv5LSJ4Q_-3nRYSqTKFyFeSgOhAL5t7pRSON0Uf5MNwuBTzrbWK46B22049dpuGS_MFhnV8cq7I-rJ-O1E4-5urvpnQ2N3kBiGoYyu5CJLy4w4Z9sxQe3SZ2_XSKga8AGmuazntpCDyDW9YFceh_h_HgpaTZtvmRX-A5ozjDdRP0NTaHY2yZTdeb_ib1H1AxKFikOHe-wWBkAFikR52VfrlCBVCK-Ni34LfYxmNxHI4qzyWIkvEssKE548zBnjVjVj3DQmfP5jTDLrFjX6FgKuUpiWF9B9ZZO-3rozcI5fF8F8ln-brt6a0KbK5Y496h9G1bD4YS6d3QhK_w9qtsTlGigbcePO7l8T3wnLgjgdDF2DpcBSbVw9h4NQs4OH1544EP7kS5yOvvsoKGJeG-H9gcEk6IlpTtaWTmYJz0TLaT6rfk3NgzF8VdWDFH9zS7CWhACIl0hPiJP1TT-x5Vb95-qWJ1H51-T6yahFyvxqxA-tDu1-ScQHrL10bE0zXuJ-ic3GqS1U-GCS6v4XxqNiNvmH_kzLwA0esZx0C8n6zzwhypcHHmjN-n6JHiPwdIJT8lCTT6KPHIv-YP_QT3tfNe2s2kTLfc25ED8PWQnYj2eshXGdiXhSZGmsgTG97lL5_Qu_s1XfYYYE5ujlCf1CYaHZBIakVRytPhpJ5JPgCb1tjm4c59PNGWwlWYXtmXqhsqCgQIaV41gMB2UZE29ZX5JK0v24qBuSTzypZ0zoey74YwCwYdrOwFI6ctPqPNk4rWhtyMFBe4dS3T5WncD03D_gIsHscwbcRWSN1ULGtLi_Ac2QagVSq24hZKnIg
dt = .1;           % sampling time
err_dist = speed_max * dt;

%%
t = 100;
t0 = t-1;

C1_old = C1(t0,:);
norm_old = norm(t0,:);
theta_old = theta(t0);
phi_old = phi(t0);

V1 = data.meas_vind_1(t);
V2 = data.meas_vind_2(t);

%%
D = sqrt(sum(C1(t,:).^2));
cond1 = abs(DIST-D) < err;

cond2 = ...
    abs(X-C1_old(1))<err_dist & ...
    abs(Y-C1_old(2))<err_dist & ...
    abs(Z-C1_old(3))<err_dist;

[~,theta_idx] = min(abs(theta_range-theta_old));
cond3 = THETA==theta_range(theta_idx) | ...
    THETA==theta_range( min(theta_idx+1,length(theta_range)) ) | ...
    THETA==theta_range( max(theta_idx-1,1) );

[~,phi_idx] = min(abs(phi_range-phi_old));
cond4 = PHI==phi_range(phi_idx) | ...
    PHI==phi_range( min(phi_idx+1,length(phi_range)) ) | ...
    PHI==phi_range( max(phi_idx-1,1) );
if phi_idx == length(phi_range)
    cond4 = cond4 + PHI==phi_range(1);
elseif phi_idx == 1
    cond4 = cond4 + PHI==phi_range(end);
end

cond5 = vind~=0 & abs(vind-V1)< err;

idx_list = find(cond1 .* cond2 .* cond3 .* cond4 .* cond5 );
sol = [];
for i = 1:length(idx_list)
    [ix, iy, iz, it, ip] = ind2sub(size(vind), idx_list(i));
    x = x_range(ix); y = y_range(iy); z = z_range(iz);
    theta = theta_range(it);  phi = phi_range(ip);
    
    n = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    c1 = [x, y, z] ;
    c2 = c1 + 3 * n;
    v2 = synthesize_vind(c2* 1e-2, n);
    
    if abs(V2 - v2) < err
        sol = [sol; x, y, z, theta, phi];
    end
end

%%
figure, 
scatter3(sol(:,1), sol(:,2), sol(:,3), 'b.')
hold on, 
scatter3(C1(t,1), C1(t,2), C1(t,3), 'ro')
xlim([x_range(1),x_range(end)]);
ylim([y_range(1),y_range(end)]);
zlim([z_range(1),z_range(end)]);

