close all; clear all; clc; addpath('src');
%% Estimate M from Vind
dataset_name = 'synth_2coils';
folder_path = strcat('../dataset/', dataset_name,'/data/');
file_list = dir(strcat(folder_path,'*.csv'));

%%
f=13.56e6;  w = 2*pi*f;
Z01=50;   Z02=1e6;

reader_coil = COIL( 15e-2, 1, 1);
reader_abcd = RESISTOR('S', 50).ABCD(f);
Ar=reader_abcd(1,1); Br=reader_abcd(1,2); 
Cr=reader_abcd(2,1); Dr=reader_abcd(2,2);
          
tag_coil = COIL( 5e-2, 1, 10);
tag_abcd = RESISTOR('S', 1).ABCD(f) * ...
    CAPACITOR('P', 1e-9).ABCD(f) * ...
    RESISTOR('P', 1e3).ABCD(f);
At=tag_abcd(1,1); Bt=tag_abcd(1,2); 
Ct=tag_abcd(2,1); Dt=tag_abcd(2,2);

a_mtx = 1/w * [Br*At, Br*Bt; Dr*At, Dr*Bt];
b_mtx = w * [Ar*Ct, Ar*Dt; Cr*Ct, Cr*Dt];
scale_mtx = [Z02, 1; Z01*Z02, Z01];
c =  2*sqrt( real(Z01) * real(Z02) ) ;
a = [1,1]*(scale_mtx.*a_mtx)*[1;1]  * 1/c;
b = [1,1]*(scale_mtx.*b_mtx)*[1;1] * 1/c;

beta = abs(a)^2 / abs(b)^2;
nwtn_step = @(m, alpha) (m^4 + alpha*m^2 + beta)/(4*m^3 + 2*alpha*m +eps); 
 
%% Synthesize
% for n = 1:length(file_list)
for n = [10]
    
    file_path = strcat(folder_path, file_list(n).name);
    data = readtable(file_path);
    Nt =  height(data);
    names = data.Properties.VariableNames;
    Ncoils = length(names(contains(names,'center')));

    M_synth = zeros(Nt,Ncoils);
    M_est = zeros(Nt,Ncoils);
    Mn1 = 1e-5;

    for t = 1:Nt
        norm = str2num(erase(data.norm{t}, ["[","]"]));        
        for i = 1:Ncoils
            center_i = str2num(erase(data.(strcat('center_',num2str(i))){t}, ["[","]"]));            
            
            if any(isnan(norm)) || any(isnan(center_i))
                M_synth(t, i) = nan;
                continue
            end
            
            tag_coil.move(center_i, norm)
            M_synth(t,i) = reader_coil.mutualInductance(tag_coil, f);
            
            vind = data.(strcat('synth_vind_',num2str(i)))(t); 
            v = (vind + .3)/5e4;
            alpha = -(a*conj(b)+conj(a)*b+1/v^2)/abs(b)^2;
    
            while true
                Mn0 = Mn1;
                Mn1 = Mn0 - nwtn_step(Mn0, alpha);               
                if abs(Mn1 - Mn0)<eps, break; end
            end
            M_est(t,i) = abs(Mn1);
        end
    end
end

%%
figure
for i = 1:Ncoils
    subplot(Ncoils,1,i), hold on
    plot(M_synth(:,i))
    plot(M_est(:,i))
end
