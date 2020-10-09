close all; clear all; clc; addpath('src');
% Generate Synthetic Vind using generated synthetic motion data
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
folder_path = strcat('../datasets/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));

%% Synthesize
N = 2;
for n = 1:length(file_list)
% for n = [1]
    file_path = strcat(folder_path, file_list(n).name);
    data = readtable(file_path);
    
    for i = 1:N
        s = data.(strcat('synth_vind_',num2str(i)));
        m = data.(strcat('meas_vind_',num2str(i)));
        
%         s = (m - mean(m))/std(m) *std(s) + mean(s);
        s(s<0) = 0; s(s>5) = 5;
        m(m<0) = 0; m(m>5) = 5;
%         
        lag = getlag(m,s);
        if lag>0
            data = data(1:end-lag,:);
            data.(strcat('meas_vind_',num2str(i))) = m(1+lag:end);
            data.(strcat('synth_vind_',num2str(i))) = s(1:end-lag);
        else
            data = data(1-lag:end,:);
            data.(strcat('meas_vind_',num2str(i))) = m(1:end+lag);
            data.(strcat('synth_vind_',num2str(i))) = s(1-lag:end);
        end
        
    end
        
    % Plot
    LW = 2;
    figure('Name','Simulation vs Measurement','NumberTitle','off');
    set(gcf, 'Units', 'Inches', 'Position', [4,1,6,4]);
    for i=1:N
        subplot(N,1,i), hold on
        plot(data.time, data.(strcat('meas_vind_',num2str(i))), ...
            'DisplayName','Measured','LineWidth', LW);
        plot(data.time, data.(strcat('synth_vind_',num2str(i))), ...
            '-.', 'DisplayName','Synthetic', 'LineWidth', LW);
        ylabel('V_{ind} (v)');
        legend('show');
        legend boxoff;
        title(strcat('Sensor_',num2str(i)));
        if i==N, xlabel('Time (s)'); end
    end
        writetable(data, file_path)
    %     saveas(gcf,strcat('figs/',num2str(n),'.png'));
end