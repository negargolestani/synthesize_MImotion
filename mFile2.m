close all; clear all; clc; addpath('src');
% Generate Synthetic Vind using generated synthetic motion data
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
% dataset_name = 'arduino';
folder_path = strcat('../dataset/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));

%% Synthesize
N = 2;
% for n = 1:length(file_list)
for n = [12]
    file_path = strcat(folder_path, file_list(n).name);
    data = readtable(file_path);
    
    meas_vind=[]; synth_vind=[]; lags=[];
    for i = 1:N
        synth_vind(i,:) = data.(strcat('synth_vind_',num2str(i)));
        synth_vind(i,synth_vind(i,:)<0) = 0;  
        synth_vind(i,synth_vind(i,:)>5) = 5;
        
        meas_vind(i,:) = data.(strcat('meas_vind_',num2str(i)));
        meas_vind(i,:) = mean(synth_vind(i,:))/mean(meas_vind(i,:)) * meas_vind(i,:);
        meas_vind(i,meas_vind(i,:)<0) = 0;  
        meas_vind(i,meas_vind(i,:)>5) = 5;  
        
        lags(i) = getlag(meas_vind(i,:), synth_vind(i,:));
        
        m = meas_vind(i,:);
        s = meas_vind(i,:);
        meas_vind(i,:) = (m - mean(m))/std(m) *std(s) + mean(s);
    end
    
    lag = round(mean(lags));

    if lag>0
        data = data(1:end-lag,:);
        for i=1:N
            data.(strcat('meas_vind_',num2str(i))) = meas_vind(i,1+lag:end)';
            data.(strcat('synth_vind_',num2str(i))) = synth_vind(i,1:end-lag)';
        end
    else
        data = data(1-lag:end,:);
        for i=1:N
            data.(strcat('meas_vind_',num2str(i))) = meas_vind(i,1:end+lag)';
            data.(strcat('synth_vind_',num2str(i))) = synth_vind(i,1-lag:end)';            
        end
    end
    
  
    % Plot
    figure
    for i=1:N
        subplot(N,1,i), hold on       
        plot(data.time, data.(strcat('meas_vind_',num2str(i))), 'DisplayName','Measured');
        plot(data.time, data.(strcat('synth_vind_',num2str(i))), 'DisplayName','Synthetic');
    ylabel('V_{ind} (V)');
    legend('show');
    title(strcat('Tag_',num2str(i)));
    if i==N, xlabel('Time (S)'); end
    end
%     writetable(data, file_path)
end