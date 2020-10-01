function [x_std] = standardize(x)
x_std = (x - mean(x)) / (std(x) + 1e-12);
end

