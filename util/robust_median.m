function med = robust_median(data, thresh, max_iter)
% function to robustly estimate the median of data, ignoring outliers
% outliers are defined at each step as data outside of the median +/- thresh*noise_std
% this is repeated for max_iter iterations

if nargin < 2, thresh = 3; end
if nargin < 3, max_iter = 5; end

mask = true(size(data)); % initialize mask to include all data
for i = 1:max_iter
    noise_std = std(data(mask)); % compute std of data, excluding outliers
    med = median(data(mask)); % compute median of data, excluding outliers
    mask = data < (med + thresh*noise_std); % mask to exclude outliers
end

end
