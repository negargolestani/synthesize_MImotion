close all; clear all; clc; addpath('src');
% Generate Synthetic Vind using motion data
%% Load motion data
% 'arduino_parallel', 'synth_parallel'
% 'arduino_orthogonal', 'synth_orthogonal'
dataset_name = 'synth_parallel';
folder_path = strcat('../datasets/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));

%% MI system
f=13.56e6;    Z01=50;   Z02=1e6;

reader_coil = COIL( 15e-2, 1, 1);
reader_abcd = RESISTOR('S', 50).ABCD(f);

tag_coil = COIL( 5e-2, 1, 10);
tag_abcd = RESISTOR('S', 1).ABCD(f) * ...
    CAPACITOR('P', 1e-9).ABCD(f) * ...
    RESISTOR('P', 1e3).ABCD(f);

%% Synthesize
for n = 1:length(file_list)
    
    file_path = strcat(folder_path, file_list(n).name);
    data = readtable(file_path);
    Nt =  height(data);
    names = data.Properties.VariableNames;
    Ncoils = sum(cell2mat(strfind(names, 'center_')));    
    synth_vind = zeros(Nt,Ncoils);
    
    for t = 1:Nt      
        for i = 1:Ncoils
            center_i = str2num(erase(data.(strcat('center_',num2str(i))){t}, ["[","]"]));
            norm_i = str2num(erase(data.(strcat('norm_',num2str(i))){t}, ["[","]"]));
            
            if any(isnan(norm_i)) || any(isnan(center_i))
                synth_vind(t, i) = nan;
                continue
            end
            tag_coil.move(center_i, norm_i)
            M = reader_coil.mutualInductance(tag_coil, f);
            ms_abcd = INDUCTOR('S',-M).ABCD(f);
            mp_abcd = INDUCTOR('P', M).ABCD(f);
            abcd = reader_abcd * ms_abcd * mp_abcd * ms_abcd * tag_abcd;
            A = abcd(1,1); B = abcd(1,2); C = abcd(2,1); D = abcd(2,2);
            S21 = (2*sqrt(real(Z01)*real(Z02)))./(A*Z02+B+C*Z01*Z02+D*Z01);
            vind = abs(S21) * 5e4 - .3;
            data.(strcat('synth_vind_',num2str(i)))(t) = vind;
        end
    end
    
    % Save
    writetable(data, file_path), n
end
