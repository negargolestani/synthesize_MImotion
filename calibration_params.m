close all; clear all; clc; addpath('src');
%% Meas vs Synth
dataset_name_list = {'arduino', 'arduino_00', 'arduino_01', ...
    'arduino_parallel', 'arduino_orthogonal'};
scale = [];
bias = [];
for j = 1:length(dataset_name_list)
    dataset_name = dataset_name_list{j};
    folder_path = strcat('../datasets/', dataset_name,'/data/');
    file_list = dir(strcat(folder_path,'*.csv'));
    
    for n = 1:length(file_list)
        file_path = strcat(folder_path, file_list(n).name);
        data = readtable(file_path);
        
        for i=1:2
            s = data.(strcat('synth_vind_',num2str(i)));
            m = data.(strcat('meas_vind_',num2str(i)));
            
            mean_s = nanmean(s); std_s = nanstd(s);
            mean_m = nanmean(m); std_m = nanstd(m); 
            
            scale(end+1) = std_s/std_m;
            bias(end+1) = mean_s - mean_m * std_s/std_m;
                       
        end
    end
end
%%
figure, plot(scale);
figure, plot(bias);

median(scale)
median(bias)
