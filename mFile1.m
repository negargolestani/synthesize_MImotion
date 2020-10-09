close all; clear all; clc; addpath('src');
% Generate Synthetic Vind using generated synthetic motion data
%% Load synthetic motion data
dataset_name = 'arduino_parallel';
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

dist = [0, 0.03];
Ncoils = length(dist);

%% Synthesize
for n = 1:length(file_list)
    
    file_path = strcat(folder_path, file_list(n).name);
    data = readtable(file_path);
    Nt =  height(data);
    names = data.Properties.VariableNames;       
%     Ncoils = length(names(contains(names,'center')));
    synth_vind = zeros(Nt,Ncoils);
    
    for t = 1:Nt
        norm = str2num(erase(data.norm{t}, ["[","]"]));
        
        for i = 1:Ncoils
            center_i = str2num(erase(data.center_1{t}, ["[","]"]));
            center_i = center_i + norm * dist(i); 
            if any(isnan(norm)) || any(isnan(center_i))
                synth_vind(t, i) = nan;
                continue
            end
            tag_coil.move(center_i, norm)
            M = reader_coil.mutualInductance(tag_coil, f);
            ms_abcd = INDUCTOR('S',-M).ABCD(f);
            mp_abcd = INDUCTOR('P', M).ABCD(f);
            abcd = reader_abcd * ms_abcd * mp_abcd * ms_abcd * tag_abcd;
            A = abcd(1,1); B = abcd(1,2); C = abcd(2,1); D = abcd(2,2);
            S21 = (2*sqrt(real(Z01)*real(Z02)))./(A*Z02+B+C*Z01*Z02+D*Z01);
            vind = abs(S21) * 5e4 - .3;
%             vind = abs(S21);
            data.(strcat('synth_vind_',num2str(i)))(t) = vind;
        end
    end
       
    % Save
    writetable(data, file_path)
    n
end
