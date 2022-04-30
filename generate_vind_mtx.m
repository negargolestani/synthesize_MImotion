close all; clear all; clc; addpath('src');
%% Generate Vind for whole space 
f=13.56e6;    Z01=50;   Z02=1e6;

reader_coil = COIL( 15e-2, 1, 1);
reader_abcd = RESISTOR('S', 50).ABCD(f);

tag_coil = COIL( 5e-2, 1, 10);
tag_abcd = RESISTOR('S', 1).ABCD(f) * ...
    CAPACITOR('P', 1e-9).ABCD(f) * ...
    RESISTOR('P', 1e3).ABCD(f);

%%
x_range = -20:20;
y_range = -20:20;
z_range = -50:-10;
theta_range = (0:5:60) * pi/180;
phi_range = (0:5:360) * pi/180;

%%
load('vind.mat');
vind_old = vind;

%%
for i = 1:length(x_range)
    for j = 1:length(y_range)
        for k = 1:length(z_range)
            for m = 1:length(theta_range)
                for n = 1:length(phi_range)
                    theta = theta_range(m);
                    phi = phi_range(n);
                    if mod(m,2)==1 && mod(n,2)==1
                        mm = (m+1)/2;   
                        nn = (n+1)/2;
                        vind(i,j,k,m,n) = vind_old(i,j,k,mm,nn);
                        continue
                    end
                    norm = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
                    center = [x_range(i), y_range(j), z_range(k)] * 1e-2;
                    tag_coil.move(center, norm);
                    M = reader_coil.mutualInductance(tag_coil, f);
                    ms_abcd = INDUCTOR('S',-M).ABCD(f);
                    mp_abcd = INDUCTOR('P', M).ABCD(f);
                    abcd = reader_abcd * ms_abcd * mp_abcd * ms_abcd * tag_abcd;
                    A = abcd(1,1); B = abcd(1,2); C = abcd(2,1); D = abcd(2,2);
                    S21 = (2*sqrt(real(Z01)*real(Z02)))./(A*Z02+B+C*Z01*Z02+D*Z01);
                    vind(i,j,k,m,n) = abs(S21) * 5e4 - .3;
                end
            end
        end
    end
    i
end
vind(vind<0) = 0;
%%

save('x_range.mat', 'x_range');
save('y_range.mat', 'y_range');
save('z_range.mat', 'z_range');
save('theta_range.mat', 'theta_range');
save('phi_range.mat', 'phi_range');
save('vind.mat', 'vind');
