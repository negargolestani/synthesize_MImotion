close all; clear all; clc; addpath('src');
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
folder_path = strcat('../dataset/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));
file_name = file_list( 1 ).name;

data = readtable( strcat(folder_path, file_name) );
data = data(50:end,:);
Nt =  height(data);
for t = 1:Nt
    C1(t,:) = str2num(erase(data.center_1{t}, ["[","]"]))*100;
    norm(t,:) = str2num(erase(data.norm{t}, ["[","]"]));
end

%%
win_Size = 5;
t = 10; 


err = 1;
speed_max = 200;    % maximum velocity (cm/s)                               ref: https://watermark.silverchair.com/M584.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAr0wggK5BgkqhkiG9w0BBwagggKqMIICpgIBADCCAp8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMnRZ6g5KHto7PK_-kAgEQgIICcOeLv5LSJ4Q_-3nRYSqTKFyFeSgOhAL5t7pRSON0Uf5MNwuBTzrbWK46B22049dpuGS_MFhnV8cq7I-rJ-O1E4-5urvpnQ2N3kBiGoYyu5CJLy4w4Z9sxQe3SZ2_XSKga8AGmuazntpCDyDW9YFceh_h_HgpaTZtvmRX-A5ozjDdRP0NTaHY2yZTdeb_ib1H1AxKFikOHe-wWBkAFikR52VfrlCBVCK-Ni34LfYxmNxHI4qzyWIkvEssKE548zBnjVjVj3DQmfP5jTDLrFjX6FgKuUpiWF9B9ZZO-3rozcI5fF8F8ln-brt6a0KbK5Y496h9G1bD4YS6d3QhK_w9qtsTlGigbcePO7l8T3wnLgjgdDF2DpcBSbVw9h4NQs4OH1544EP7kS5yOvvsoKGJeG-H9gcEk6IlpTtaWTmYJz0TLaT6rfk3NgzF8VdWDFH9zS7CWhACIl0hPiJP1TT-x5Vb95-qWJ1H51-T6yahFyvxqxA-tDu1-ScQHrL10bE0zXuJ-ic3GqS1U-GCS6v4XxqNiNvmH_kzLwA0esZx0C8n6zzwhypcHHmjN-n6JHiPwdIJT8lCTT6KPHIv-YP_QT3tfNe2s2kTLfc25ED8PWQnYj2eshXGdiXhSZGmsgTG97lL5_Qu_s1XfYYYE5ujlCf1CYaHZBIakVRytPhpJ5JPgCb1tjm4c59PNGWwlWYXtmXqhsqCgQIaV41gMB2UZE29ZX5JK0v24qBuSTzypZ0zoey74YwCwYdrOwFI6ctPqPNk4rWhtyMFBe4dS3T5WncD03D_gIsHscwbcRWSN1ULGtLi_Ac2QagVSq24hZKnIg
Ts = .01;           % sampling time
er_dist = speed_max * Ts;
% win_Size = size(C1_0,1);

load('vind.mat');
load('x_range.mat');
load('y_range.mat');
load('z_range.mat');
load('theta_range.mat');
load('phi_range.mat');
[X, Y, Z, THETA, PHI] = ndgrid(x_range, y_range, z_range, theta_range, phi_range);

% Inputs
V1 = data.meas_vind_1(t);
V2 = data.meas_vind_2(t);

C1_old = median(C1(t-win_Size:t-1,:), 1);
norm_old = median(norm(t-win_Size:t-1,:), 1);

%%
dist = sqrt( (X-C1_old(1)).^2 + (Y-C1_old(2)).^2 + (Z-C1_old(3)).^2 );
cond1 = dist < er_dist;

theta = acos(norm_old(3));
[~,theta_idx] = min(abs(theta_range-theta));
cond2 = THETA==theta_range(theta_idx) | ...
    THETA==theta_range( min(theta_idx+1,length(theta_range)) ) | ...
    THETA==theta_range( max(theta_idx-1,1) );

phi = mod(atan2(norm_old(2),norm_old(1)), 2*pi);
[~,phi_idx] = min(abs(phi_range-phi));
cond3 = PHI==phi_range(phi_idx) | ...
    PHI==phi_range( min(phi_idx+1,length(phi_range)) ) | ...
    PHI==phi_range( max(phi_idx-1,1) );
if phi_idx == length(phi_range)
    cond3 = cond3 + PHI==phi_range(1);
elseif phi_idx == 1
    cond3 = cond3 + PHI==phi_range(end);
end

cond4 = vind~=0 & (abs(vind-V1)< err | abs(gradient(sign(vind-V1)))==2);

dist = sqrt( X.^2 + Y.^2 + Z.^2 );
cond5 = sqrt(sum(C1(t,:).^2))- 1 < dist &  dist < sqrt(sum(C1(t,:).^2))+1;

cond6 = abs(gradient(vind)) < abs(V2-V1) + .3;


%%
idx_list = find(cond1 .* cond2 .* cond3 .* cond4 .* cond5 .*cond6);

figure,
for i = 1:length(idx_list)
    idx = idx_list(i);
%     
    [ix, iy, iz, it, ip] = ind2sub(size(vind), idx);
    x = x_range(ix); 
    y = y_range(iy); z = z_range(iz);
%     theta = theta_range(it);  phi = phi_range(ip);
    scatter3(x,y,z,'.'); hold on
end
scatter3(C1(t,1), C1(t,2), C1(t,3),'ro');
xlim([-20, 20]); 
ylim([-20, 20]);
zlim([-50, -10]);