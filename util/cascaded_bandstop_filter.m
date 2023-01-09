function y = cascaded_bandstop_filter(x, fs, f0, N, bw0, fmin, fmax)
    % x = input data
    % N = filter order
    % f0 = fundamental frequency to filter (e.g. 60 Hz)
    % fs = sampling rate
    % bw_c = bandwidth of fundamental bandstop filter
    % fmin = minimum frequency to filter
    % fmax = maximum frequency to filter
    
    if nargin < 3, f0 = 60; end
    if nargin < 4, N = 20; end
    if nargin < 5, bw0 = 1; end
    if nargin < 6, fmin = 0; end
    if nargin < 7, fmax = fs/2 - bw0/2*(fs/2)/f0; end

    freqs = f0:f0:fs/2;
    freqs(freqs < fmin | freqs > fmax) = [];
    
    Hd_ca = [];
    eval_str = 'dfilt.cascade(';
    for i = 1:length(freqs)
        freq = freqs(i);
%         fc1 = freq - bw0/2*freq/f0;  % First Cutoff Frequency
%         fc2 = freq + bw0/2*freq/f0;  % Second Cutoff Frequency
        fc1 = freq - bw0/2;  % First Cutoff Frequency
        fc2 = freq + bw0/2;  % Second Cutoff Frequency
    
        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandstop('N,F3dB1,F3dB2', N, fc1, fc2, fs);
        Hd = design(h, 'butter');
        Hd_ca{i} = Hd;
        eval_str = [eval_str, sprintf('Hd_ca{%d},', i)];
    end
    eval_str = [eval_str(1:end-1), ')'];
    
    Hd = eval(eval_str);

    y = filter(Hd, x);
end