function [C1_pred_, norm_pred_] = estimate(C1_0, norm_0, V1_list, V2_list)

err = 1;
speed_max = 200;    % maximum velocity (cm/s)                               ref: https://watermark.silverchair.com/M584.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAr0wggK5BgkqhkiG9w0BBwagggKqMIICpgIBADCCAp8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMnRZ6g5KHto7PK_-kAgEQgIICcOeLv5LSJ4Q_-3nRYSqTKFyFeSgOhAL5t7pRSON0Uf5MNwuBTzrbWK46B22049dpuGS_MFhnV8cq7I-rJ-O1E4-5urvpnQ2N3kBiGoYyu5CJLy4w4Z9sxQe3SZ2_XSKga8AGmuazntpCDyDW9YFceh_h_HgpaTZtvmRX-A5ozjDdRP0NTaHY2yZTdeb_ib1H1AxKFikOHe-wWBkAFikR52VfrlCBVCK-Ni34LfYxmNxHI4qzyWIkvEssKE548zBnjVjVj3DQmfP5jTDLrFjX6FgKuUpiWF9B9ZZO-3rozcI5fF8F8ln-brt6a0KbK5Y496h9G1bD4YS6d3QhK_w9qtsTlGigbcePO7l8T3wnLgjgdDF2DpcBSbVw9h4NQs4OH1544EP7kS5yOvvsoKGJeG-H9gcEk6IlpTtaWTmYJz0TLaT6rfk3NgzF8VdWDFH9zS7CWhACIl0hPiJP1TT-x5Vb95-qWJ1H51-T6yahFyvxqxA-tDu1-ScQHrL10bE0zXuJ-ic3GqS1U-GCS6v4XxqNiNvmH_kzLwA0esZx0C8n6zzwhypcHHmjN-n6JHiPwdIJT8lCTT6KPHIv-YP_QT3tfNe2s2kTLfc25ED8PWQnYj2eshXGdiXhSZGmsgTG97lL5_Qu_s1XfYYYE5ujlCf1CYaHZBIakVRytPhpJ5JPgCb1tjm4c59PNGWwlWYXtmXqhsqCgQIaV41gMB2UZE29ZX5JK0v24qBuSTzypZ0zoey74YwCwYdrOwFI6ctPqPNk4rWhtyMFBe4dS3T5WncD03D_gIsHscwbcRWSN1ULGtLi_Ac2QagVSq24hZKnIg
Ts = .01;           % sampling time
er_dist = speed_max * Ts;
% win_Size = size(C1_0,1);
win_Size = 5;

load('vind.mat');
load('x_range.mat');
load('y_range.mat');
load('z_range.mat');
load('theta_range.mat');
load('phi_range.mat');
[X, Y, Z, THETA, PHI] = ndgrid(x_range, y_range, z_range, theta_range, phi_range);

Nt = length(V1_list);
C1_pred(1:win_Size,:) = C1_0(1:win_Size,:);
norm_pred(1:win_Size,:) = norm_0(1:win_Size,:);

for t = win_Size+1:Nt
    V1 = V1_list(t);
    V2 = V2_list(t);
    
    C1_old = median(C1_pred(t-win_Size:t-1,:), 1);
    norm_old = median(norm_pred(t-win_Size:t-1,:), 1);
    
    dist = sqrt( (X-C1_old(1)).^2 + (Y-C1_old(2)).^2 + (Z-C1_old(3)).^2 );
    cond1 = dist < er_dist;
    
    dist = sqrt( X.^2 + Y.^2 + Z.^2 );
    cond5 = sqrt(sum(C1_0(t,:).^2))-2 < dist &  dist < sqrt(sum(C1_0(t,:).^2))+2;
    
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
    idx_list = find(cond1 .* cond2 .* cond3 .* cond4 .* cond5);
    
    
    C1_pred_list = [];
    norm_pred_list = [];
    for i = 1:length(idx_list)
        idx = idx_list(i);
        
        [ix, iy, iz, it, ip] = ind2sub(size(vind), idx);
        x = x_range(ix); y = y_range(iy); z = z_range(iz);
        theta = theta_range(it);  phi = phi_range(ip);
        
        norm = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        c1 = [x,y,z] ;
        c2 = get_center2(c1,norm);
        
        v2_f = vind( ...
            find(x_range==floor(c2(1))),...
            find(y_range==floor(c2(2))),...
            find(z_range==floor(c2(3))),...
            it,ip);
        v2_c = vind( ...
            find(x_range==ceil(c2(1))),...
            find(y_range==ceil(c2(2))),...
            find(z_range==ceil(c2(3))),...
            it,ip);
        %         if min(v2_f, v2_c)<V2 & V2<max(v2_f,v2_c)
        %             v2 = synthesize_vind(c2*1e-2, norm);
        %             if abs(v2 - V2) > err
        %                 continue;
        %             end
        %         else
        %             continue;
        %         end
        if  V2<min(v2_f,v2_c)-err | V2>max(v2_f,v2_c)+err
            continue;
        end
        %     vp = get_Vdphi(c1,c2);
        %     if abs(vp - Vp) > err
        %         continue;
        %     end
        
        C1_pred_list = [C1_pred_list; c1];
        norm_pred_list = [norm_pred_list; norm];
    end
    
    if any(C1_pred_list)
        C1_pred(t,:) = median(C1_pred_list, 1);
    else
        C1_pred(t,:) = C1_old;
    end
    if any(norm_pred_list)
        norm_pred(t,:) = median(norm_pred_list, 1);
    else
        norm_pred(t,:) = norm_old;
    end
    
    C1_pred_{t} = C1_pred_list;
    norm_pred_{t} = norm_pred_list;
    t
end

end