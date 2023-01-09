function baseline = smooth_masked_fit(data, mask, lambda)
% estimate baseline as a smoothed fit to masked data (outliers excluded)
% data is the input data to fit
% mask is the logical mask of valid indices
% lambda is the smoothing parameter

if nargin < 3, lambda = 1e12; end % default lambda (works well for example data)
if nargin < 2, mask = true(size(data)); end % default to use all data

N = length(data);
D = diff(speye(N), 2);
w = double(mask);

W = spdiags(w, 0, N, N);
C = chol(W + lambda * (D.' * D));
baseline = C \ (C.' \ (w .* data));

% figure(99); clf; hold on;
% temp = data;
% temp(~mask) = NaN;
% plot(data);
% plot(temp, 'linewidth', 5);
% plot(baseline);

end