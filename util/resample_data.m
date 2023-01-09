function data_resampled = resample_data(data, samplerate, resample_rate)
% resamples data with padding to reduce edge effects
% Inputs:
% data = nx1 column vector of data values (any units)
% samplerate = sample rate of data, in Hz (must be integer)
% resample_rate = resample rate of data, in Hz (THIS SHOULD BE AN INTEGER TO AVOID SUS BEHAVIOR)
% Outputs:
% data_resampled = mx1 column vector of resampled data values (same units)

if mod(resample_rate,1)~=0
    error('resample_rate must be an integer! asked for %f', resample_rate);
end

[data_pad, NR] = pad_data(data, samplerate); % pad data to reduce edge effects

% % might speed up resampling
% % requires samplerate & resample_rate to be integers
% g = gcd(samplerate, resample_rate);
% data_resampled_pad = resample(data_pad, resample_rate/g, samplerate/g);

data_resampled_pad = resample(data_pad, resample_rate, samplerate);

nr = floor( NR * resample_rate / samplerate );
m = ceil(length(data) / samplerate * resample_rate); % length of output of resample.m
data_resampled = data_resampled_pad((1:m)+nr);

end