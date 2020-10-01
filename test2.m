close all; clear all; clc; addpath('src');
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
folder_path = strcat('../dataset/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));
file_name = file_list( 1 ).name;

data = readtable( strcat(folder_path, file_name) );
data = data(50:end,:);
Nt = height(data);
for t = 1:Nt
    C1(t,:) = str2num(erase(data.center_1{t}, ["[","]"]))*100;
    norm(t,:) = str2num(erase(data.norm{t}, ["[","]"]));
end

[C1_pred_, norm_pred_] = estimate( ...
    C1, ...
    norm, ...
    data.synth_vind_1, ...
    data.synth_vind_2 ...
    );
%% Plot
figure, set(gcf, 'Units', 'Inches', 'Position', [3,1,9,7]);
for i = 1:3
    subplot(3,1,i);    hold on
    plot(C1(:,i), 'Displayname', 'Target')
    for t = 1:Nt
        c1p = C1_pred_{t};
        if any(c1p)
            c1p = c1p(:,i);
            scatter(t*ones(size(c1p)), c1p, 'r.');
        end
    end
end

%% Plot
% C1_pred_ = zeros(size(C1));
% for t = 1:length(C1_pred)
%     if any(C1_pred{t})
%         C1_pred_(t,:) = median(C1_pred{t},1);
%     end
% end
% % C1_pred_ = smoothdata(C1_pred_);
% figure, set(gcf, 'Units', 'Inches', 'Position', [3,1,9,7]);
% for i = 1:3
%     subplot(3,1,i);    hold on
%     plot(C1(:,i), 'Displayname', 'Target')
%     C1_pred_(:,i) = smoothdata(C1_pred_(:,i));
%     plot(C1_pred_(:,i), 'Displayname', 'Estimation')    
% end