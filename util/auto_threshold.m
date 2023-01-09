function thresh = auto_threshold(data, x0, x1)
% generates a guess for a threshold based on finding knee of qqplot
% assumes qqplot looks like __._/, where . is the origin
% x0 is the starting standard quantile for the kneedle algorithm
% x1 is the ending standard quantile for the kneedle algorithm

if nargin < 3, x1 = 2; end % default q1 of 2 std above mean (captures 95% of data for normal)
if nargin < 2, x0 = -2; end % default q0 of 2 std below mean (captures 95% of data for normal)

qq = qqplot(data); % compute quantile-quantile plot
x = qq(1).XData.'; % x coords of qq plot
y = qq(1).YData.'; % y coords of qq plot
close(gcf);

mask = x >= x0 & x <= x1; % mask of valid indices to use for kneedle algorithm
x = x(mask);
y = y(mask);

y01 = interp1(x, y, [x0, x1], [], 'extrap');
y0 = y01(1);
y1 = y01(2);

kneedle = abs((y1-y0).*x - (x1-x0).*y + x1*y0 - y1*x0) ./ sqrt((y1-y0)^2 + (x1-x0)^2); 

[~, ind] = max(kneedle);
if ind == length(kneedle)
    ind = ind - 1; 
end

thresh = y(ind);

end