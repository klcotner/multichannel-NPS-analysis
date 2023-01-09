function [norm_corr_map, corr_map] = apply_FB(input, FB, use_unit_norm)
%{

INPUTS
* input = (1xn row vector of data)
    > this is usually a block of data to be filtered
    > n should be greater than m (the filter length)
* FB = (num_filts x m filterbank)
    > each row is a different filter, but all are length m
    > m should be less than n (the input length)
* use_unit_norm = (logical flag for using unit norm filters)
    > if true, this function normalizes each filter to unit norm before
      performing cross-correlation
    > defaults to true

OUTPUTS
* norm_corr_map = (num_filts x n array of normalized correlation coefficients)
* corr_map = (num_filts x n array of filtered data a.k.a. correlation map)
%}

if nargin < 3
    use_unit_norm = true; % default to unit norm filters
end
if size(input,2)==1
    if size(input,1)>1 % input is a column vector
        input = input.';
    else
        error('input must be a vector');
    end
end

[num_filts, m] = size(FB);
n = length(input);
if m >= n
    error('m (the filter length) must be less than n (the input length)');
end

norm_corr_map = nan(num_filts, n+m-1);
corr_map = nan(num_filts, n+m-1); % initialize output

for i = 1:size(FB,1)
    filt = FB(i,:);
    if use_unit_norm
        filt = filt ./ norm(FB(i,:));
    end
    [nxc, xc] = normxcorr(filt, input);
    norm_corr_map(i,:) = nxc;
    corr_map(i,:) = xc;
end

offset = ceil((m-1)/2); % offset to remove (m-1) edge values
norm_corr_map = norm_corr_map(:, offset + (1:n));
corr_map = corr_map(:, offset + (1:n)); % trim output

end

