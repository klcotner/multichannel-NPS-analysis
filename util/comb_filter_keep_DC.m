function Y = comb_filter_keep_DC(X, f0, fs, bw_c, fmin, fmax)
% this function does comb filtering without removing the DC component
% X is the input signal
% f0 is the fundamental frequency to filter out (e.g. 60Hz) [Hz]
% fs is the sampling rate of the input signal [samp/s]
% bw_c is the bandwidth of each null [Hz] (default 1)
% fmin is the minimum frequency (exclusive) [Hz] (default 0)
% fmax is the maximum frequency (inclusive) [Hz] (default fs/2)
% Y is the comb-filtered output

if nargin < 4, bw_c = 1; end
if nargin < 5, fmin = 0; fmax = fs/2; end % default to Nyquist frequency
if isempty(bw_c), bw_c = 1; end

n = fs/f0;
if mod(n,1)~=0
    error('Comb_filter_keep_DC cannot filter this frequency at this sampling rate!');
end
[b, a] = iircomb(n, bw_c/(fs/2));
[z, p, k] = tf2zpk(b, a);
freqs_z = abs(angle(z) / pi * fs/2);
mask_z = freqs_z > fmin & freqs_z <= fmax;
freqs_p = abs(angle(p) / pi * fs/2);
mask_p = freqs_p > fmin & freqs_p <= fmax;
[sos, g] = zp2sos(z(mask_z), p(mask_p), k);
Y = filtfilt(sos, g, X-X(1,:)) + X(1,:);

end

