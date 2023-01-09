function noise_std = robust_std(data, thresh, max_iter)
% function to robustly estimate the noise std of data, ignoring outliers
% outliers are defined at each step as data outside of the mean +/- thresh*noise_std
% this is repeated for max_iter iterations
% note that the baseline must already have been subtracted for this to work

if nargin < 2, thresh = 3; end
if nargin < 3, max_iter = 5; end

mask = true(size(data)); % initialize mask to include all data
for i = 1:max_iter
    noise_std = std(data(mask)); % compute std of data, excluding outliers
    mask = data < (mean(data(mask)) + thresh*noise_std); % mask to exclude outliers
end

end
