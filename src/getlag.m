function [ lag ] = getlag(a, b)

corr_max = 0;
lag = 0; 
for d = -50:50
    if d>0
        a_sh = a(1+d:end);
        b_sh = b(1:end-d);
    else
        a_sh = a(1:end+d);
        b_sh = b(1-d:end);
    end
    
    c = corrcoef(a_sh, b_sh, 'rows','complete');
    c = c(1,2);
    if c > corr_max
        corr_max = c;
        lag = d;
    end       
end

end


