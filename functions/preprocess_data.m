% INPUTS:
% * raw_data = nx1 column vector of doubles: raw resistance values [ohms]
% * raw_time = nx1 column vector of doubles: raw time values [sec]
% * sample_rate = acquisition sample rate [Hz]
% * FB_params (struct) = filter bank information (to infer pulse widths) with fields:
%     > sequence_len = int = number of segments (1's and 0's) in the encoded sequence(s)
%     > base_time (optional) = expected transit time of one event [sec] (default = 0.2 sec)
%     > min_ratio (optional) = shortest signal template to generate, as a fraction of base_time (default = 0.05)
%     > max_ratio (optional) = longest signal template to generate, as a fraction of base_time (default = 1.5)
% * pp_params (optional, struct) = preprocessing parameters structure with the following fields:
%     > resample_rate = desired resampling rate [Hz]
%         - default = 5 * bw_lpf
%         - if apply_60Hz_comb, will overwrite with next-highest multiple of 60
%     > apply_60Hz_comb = whether to apply comb filtering
%         - default = false
%     > bw_lpf = lowpass filter bandwidth (cutoff) [Hz]
%         - default = 5/min_pulse_width
%     > lambda = smoothing parameter for baseline fit
%         - default = 1e12
%     > bw_c = comb filter bandwidth [Hz]
%         - default = 1
% * plotting (optional) = bool
%	  - if true, plots the raw data, filtered data, baseline fit, and noise mask
%	  - default = false
% * new_filepath (optional) = bool or string/char
%	  - if true, saves preprocessed data to "preprocessed_data.mat" in the current directory
%	  - if string/char, saves preprocessed data to new_filepath (appends ".mat" if necessary)
%	  - if false/empty, doesn't save the preprocessed data
%     - default = false
% 
% OUTPUTS:
% * resampled_time = mx1 column vector of doubles : time values [sec] at the resampled rate
% * filtered_data = mx1 column vector of doubles : filtered resistance measurements [ohms] at the resampled rate
% * baseline_data = mx1 column vector of doubles : initial baseline fit [ohms] at the resampled rate
% * noise_mask = mx1 column vector of logicals : mask that is true where the data is only noise, not signal
% * resample_rate = double : the actual resample rate used [Hz]
% * pp_params_used = struct : the final values used for each preprocessing step with fields:
%     > resample_rate = sample rate of filtered data [Hz]
%     > apply_60Hz_comb = whether comb filtering was used
%     > bw_lpf = lowpass filter bandwidth (cutoff) [Hz]
%     > lambda = smoothing parameter for baseline fit
%     > min_pulse_width = minimum expected single-pulse width [sec]
%     > max_pulse_width = maximum expected single-pulse width [sec]
%     > med_filtlen = width of median filter used [samples]
%     > minmax_filtlen = window width for moving min&max used for envelope calculation [samples]
%     > env_thresh = difference threshold used to generate noise mask [ohms]
% 
% NOTES:
% * TODO: steepness parameter in call to lowpass
% * this function does NOT subtract the baseline from resampled_data

function [resampled_time, filtered_data, baseline_est, noise_mask, resample_rate, pp_params_used] = ...
    preprocess_data(raw_data, raw_time, sample_rate, FB_params, pp_params, plotting, new_filepath)

%% Parse input arguments

if nargin<7
    new_filepath = false;
end
if nargin<6
    plotting = false;
end
if nargin<5
    pp_params = [];
end

if isempty(plotting)
    plotting = false;
end
if isempty(new_filepath)
    new_filepath = false;
end

%% check shape of data vectors

if isvector(raw_data)
    if isrow(raw_data)
        warning('transposing raw_data to create column vector');
        raw_data = raw_data';
    end
else
    error('raw_data must be a vector!');
end

if isvector(raw_time)
    if isrow(raw_time)
        warning('transposing raw_time to create column vector');
        raw_time = raw_time';
    end
else
    error('raw_time must be a vector!');
end

if ~all(size(raw_data)==size(raw_time))
    error('raw_data and raw_time must be the same length!');
end

%% check other inputs

% sample_rate
if ~isnumeric(sample_rate) || length(sample_rate)~=1
    error('invalid input for sample_rate');
end

% FB_params
if ~isstruct(FB_params)
    error('invalid input for FB_params');
end

% pp_params
if ~isempty(pp_params) && ~isstruct(pp_params)
    error('invalid input for pp_params');
end

% plotting flag
if ~islogical(plotting)
    if ~isnumeric(plotting) || isnan(plotting)
        error('invalid input for plotting');
    end
end

% new_filepath (saving flag)
if ~islogical(new_filepath)

    % convert numeric to logical, if needed
    if isnumeric(new_filepath)
        if isnan(new_filepath)
            error('invalid input for new_filepath');
        else
            new_filepath = new_filepath~=0;
        end

    % otherwise, check for string/char
    elseif ~ischar(new_filepath) && ~isstring(new_filepath)
        error('invalid input for new_filepath');
    end

end

%% fill in default field values for FB_params

if isempty(FB_params)
    FB_params = struct();
end

if ~isfield(FB_params, 'sequence_len') || isempty(FB_params.sequence_len)
    error('sequence_len must be provided in FB_params');
elseif ~isnumeric(FB_params.sequence_len) || length(FB_params.sequence_len)~=1 || isnan(FB_params.sequence_len)
    error('invalid input for FB_params.sequence_len');
end

if ~isfield(FB_params, 'base_time') || isempty(FB_params.base_time)
    FB_params.base_time = 0.2; % [sec]
end

if ~isfield(FB_params, 'min_ratio') || isempty(FB_params.min_ratio)
    FB_params.min_ratio = 0.05;
end

if ~isfield(FB_params, 'max_ratio') || isempty(FB_params.max_ratio)
    FB_params.max_ratio = 1.5;
end

%% fill in default field values for pp_params

if isempty(pp_params)
    pp_params = struct();
end

if ~isfield(pp_params, 'resample_rate')
    pp_params.resample_rate = [];
end

if ~isfield(pp_params, 'apply_60Hz_comb') || isempty(pp_params.apply_60Hz_comb)
    pp_params.apply_60Hz_comb = false;
end

if ~isfield(pp_params, 'bw_lpf')
    pp_params.bw_lpf = [];
end

if ~isfield(pp_params, 'lambda') || isempty(pp_params.lambda)
    pp_params.lambda = 1e12;
end

if ~isfield(pp_params, 'bw_c') || isempty(pp_params.bw_c)
    pp_params.bw_c = 1;
end

%% Parse FB_params
% infer minimum & maximum expected pulse widths

base_pulse_width = FB_params.base_time / FB_params.sequence_len; % [sec]
min_pulse_width = base_pulse_width * FB_params.min_ratio; % [sec]
max_pulse_width = base_pulse_width * FB_params.max_ratio; % [sec]

%% Parse pp_params

% check resample_rate compatibility with 60Hz comb filter, if both were specified
if pp_params.apply_60Hz_comb && ~isempty(pp_params.resample_rate) && mod(pp_params.resample_rate,60)~=0
    pp_params.resample_rate = floor(pp_params.resample_rate / 60) * 60 + 60;
    warning('you specified apply_60Hz_comb as well as a resample_rate that is not a multiple of 60. Resampling rate set to %d Hz.', pp_params.resample_rate);
end

% if needed, calculate bw_lpf based on pulse_width
if isempty(pp_params.bw_lpf)
    pp_params.bw_lpf = 5/min_pulse_width;
    fprintf('Inferred bw_lpf from min_pulse_width: %f Hz\n', pp_params.bw_lpf);
end

% if needed, infer resample_rate from bw_lpf
if isempty(pp_params.resample_rate)

    % downsampled rate should be at least 5X Nyquist
    pp_params.resample_rate = 2 * pp_params.bw_lpf * 5; % [Hz]

    % if using 60hz comb filter, downsampled rate should be divisible by 60
    if pp_params.apply_60Hz_comb
        pp_params.resample_rate = floor(pp_params.resample_rate / 60) * 60 + (mod(pp_params.resample_rate, 60)~=0) * 60;
    end

    fprintf('Inferred resample_rate from bw_lpf: %d Hz\n', pp_params.resample_rate);

end

%% Resample data to lower rate (speeds up next steps)
% uses matlab resample function, which applies an antialiasing LPF

resample_rate = pp_params.resample_rate;

resampled_data = unpad_data(resample(pad_data(raw_data-raw_data(1), sample_rate), resample_rate, sample_rate), resample_rate) + raw_data(1); % [ohms]
resampled_time = (0:length(resampled_data)-1).' / resample_rate + raw_time(1); % [sec]

%% Remove 60 Hz & harmonics with comb filter
if pp_params.apply_60Hz_comb
    combed_data = unpad_data(comb_filter_keep_DC(pad_data(resampled_data-resampled_data(1), resample_rate), 60, resample_rate, pp_params.bw_c), resample_rate) + resampled_data(1);
else
    combed_data = resampled_data;
end

%% Low pass filter
filtered_data = unpad_data(lowpass(pad_data(combed_data-combed_data(1), resample_rate), pp_params.bw_lpf, resample_rate), resample_rate) + combed_data(1);

%% Median filter
% for extra smoothing without losing edges
med_filtlen = floor(max_pulse_width * resample_rate); % length of median filter [samples]
medfilted_data = unpad_data(movmedian(pad_data(filtered_data, resample_rate), med_filtlen), resample_rate);

%% Compute signal envelope

minmax_filtlen = floor(max_pulse_width * resample_rate * 4); % length of movmax and movmin filters [samples]

ub = unpad_data(movmax(pad_data(medfilted_data, resample_rate), minmax_filtlen), resample_rate); % upper bound of signal "envelope"
lb = unpad_data(movmin(pad_data(medfilted_data, resample_rate), minmax_filtlen), resample_rate); % lower bound of signal "envelope"

env = ub - lb; % signal envelope

%% Compute noise mask
env_thresh = auto_threshold(env); % compute automatic threshold to reject transit events
noise_mask = env < env_thresh; % mask of where the data is just noise (no transit events)

%% Estimate baseline using the noise mask
baseline_est = smooth_masked_fit(medfilted_data, noise_mask, pp_params.lambda);

%% save final parameters to output struct

pp_params_used = pp_params;

pp_params_used.min_pulse_width = min_pulse_width;
pp_params_used.max_pulse_width = max_pulse_width;

pp_params_used.med_filtlen = med_filtlen;
pp_params_used.minmax_filtlen = minmax_filtlen;

pp_params_used.env_thresh = env_thresh;

%% Plot (if required)

if plotting
    noise_data = medfilted_data;
    noise_data(~noise_mask) = NaN;
    
    figure;
    hold on;
    plot(raw_time, raw_data, 'DisplayName', 'Raw Data');
    plot(resampled_time, filtered_data, 'linewidth', 1.5, 'DisplayName', 'Filtered Data');
    plot(resampled_time, noise_data, 'linewidth', 2.5, 'DisplayName', 'Detected Noise');
    plot(resampled_time, baseline_est, 'linewidth', 2, 'DisplayName', 'Baseline Estimate');
    xlabel('Time [s]');
	ylabel('Amplitude [a.u.]');
	title('Preprocessed Data');
    axis tight;
    lgd = legend();
    set(lgd, 'location','southeast');

end

%% Save preprocessed data (if required)

if islogical(new_filepath) && new_filepath
    new_filepath = 'preprocessed_data.mat';
end
if ischar(new_filepath)
    if ~endsWith(new_filepath, '.mat')
        new_filepath = [new_filepath, '.mat'];
    end

    S_out = struct();
    S_out.resampled_time = resampled_time;
    S_out.filtered_data = filtered_data;
    S_out.baseline_est = baseline_est;
    S_out.noise_mask = noise_mask;
    S_out.resample_rate = resample_rate;
    S_out.pp_params_used = pp_params_used;
end

end
