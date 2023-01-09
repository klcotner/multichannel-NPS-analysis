%{
[bFB_hr, uFB_hr, index_dict_hr, bFB_lr, uFB_lr, index_dict_lr, transit_times, bw_list] = generate_FB(sequence, sample_rate, resample_rate, [FB_param], [plotting])

INPUTS
* sequence = (1xn row vector of 1's and 0's)
    > the sequence that is encoded in the device geometry (1=pore, 0=node)
* samplerate = acquisition sample rate, in Hz
* resample_rate = desired sample rate after resampling, in Hz
* FB_param = (optional, struct) = filter bank information
    > to use non-default values, may include fields:
        * num_filters = number of signal templates to generate in the filterbank (default = 100)
        * base_time = expected transit time of one event [sec] (default = 0.2 sec)
        * min_ratio = shortest signal template to generate, as a fraction of base_time (default = 0.05)
        * max_ratio = longest signal template to generate, as a fraction of base_time (default = 1.5)
	> if FB_param is not provided, default values are used for all 4 parameters
* plotting = (optional, bool) = whether or not to plot the filterbank
    > default = false

OUTPUTS
* _hr suffix = high-resolution (original) sample rate
* _lr suffix = low-resolution (downsampled) sample rate
* bFB_ = mxs array of doubles = bipolar filterbank
    > each row is a template with a different transit time
    > used for correlation; template consists of positive & negative numbers
    > magnitude is scaled such that norm(template)=1
* uFB_ = mxs array of 1's and 0's = unipolar filterbank
    > each row is a template with a different transit time
    > used for fitting; template consists of 1's and 0's
* index_dict_ = mx3 array = index information for the filterbanks
    > [length_in_samples, start_index, stop_index] for the corresponding filterbank rows
* transit_times = mx1 column vector of doubles
    > transit times of the corresponding filterbank rows (in seconds)
* bw_list = mx1 column vector of doubles
    > a rough measure of bandwith for the corresponding filterbank rows
    > units of sample^-1 (at the resampled rate)

NOTES
>> todo: support for calibrated rcode
>> note: '../Util' must be on the search path
    > for check_bool, check_sequence, and plot_FB
%}

function [bFB_hr, uFB_hr, index_dict_hr, bFB_lr, uFB_lr, index_dict_lr, transit_times, bw_list] = generate_FB(sequence, sample_rate, resample_rate, FB_param, plotting)

    %% input parser
    p = inputParser;
    p.FunctionName = 'generate_FB';
    p.PartialMatching = false;
    p.CaseSensitive = true;
    p.StructExpand = false;
    
    % required parameters
    p.addRequired('sequence', @check_sequence);
    p.addRequired('sample_rate', @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addRequired('resample_rate', @(x) isnumeric(x) && isscalar(x) && x>0);
    
    % optional parameters
    p.addOptional('FB_param', struct());
    p.addOptional('plotting', false, @check_bool);

    % run parser based on the provided arguments
    if nargin==3
        parse(p, sequence, sample_rate, resample_rate);
    elseif nargin==4
        parse(p, sequence, sample_rate, resample_rate, FB_param);
    else
        parse(p, sequence, sample_rate, resample_rate, FB_param, plotting);
    end
    
    % get plotting flag
    plotting = p.Results.plotting;

    %% parse optional structure argument

    FB_param = p.Results.FB_param;

    p2 = inputParser;
    p2.FunctionName = 'generate_FB';
    p2.PartialMatching = false;
    p2.CaseSensitive = true;
    p2.StructExpand = true;

    % optional parameters (from the FB_param struct)
    p2.addOptional('num_filters', 100, @(x) isnumeric(x) && isscalar(x) && x>0);
    p2.addOptional('min_ratio', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0);
    p2.addOptional('max_ratio', 1.5, @(x) isnumeric(x) && isscalar(x) && x>0);
    p2.addOptional('base_time', 0.2, @(x) isnumeric(x) && isscalar(x) && x>0);

    % run parser
    parse(p2, FB_param);

    % get user-defined or default values for the filter bank from FB_param
    num_filters = p2.Results.num_filters;
    base_time = p2.Results.base_time;
    min_ratio = p2.Results.min_ratio;
    max_ratio = p2.Results.max_ratio;
    
    %% construct list of transit times & base signal templates
    
    % generate list of transit times [signal lengths for the FB]
    transit_times = linspace(min_ratio*base_time, max_ratio*base_time, num_filters)'; % [sec]
    transit_times_hr = round(transit_times * sample_rate); % [# of HR samples]
    transit_times_lr = round(transit_times * resample_rate); % [# of LR samples]
    
    % maximum template length in samples at each resolution
    maxtime_hr = max(transit_times_hr); % [# of HR samples]
    maxtime_lr = max(transit_times_lr); % [# of LR samples]
    
    % find max transit time in samples where each segment has equal length
    num_seg = length(sequence); % # of segments in the sequence
    segment_length = ceil( maxtime_hr / num_seg ); % [# of samples]
    
    % construct base template sequence (unipolar)
    %   length >= the longest high-resolution transit time, in samples
    %   each segment should have the same number of samples
    basetemplate_up = kron(sequence, ones(1,segment_length));
    
    %% generate high resolution filter banks (original sample rate)
    
    % initialize filterbanks, bw_list, & index_dict (high resolution)
    % each row corresponds to a different transit time
    bFB_hr = zeros(num_filters, maxtime_hr); % each row is a template
    uFB_hr = bFB_hr; % each row is a template
    bw_list = zeros(num_filters,1); % column vector
    index_dict_hr = zeros(num_filters, 3); % [ {template length in # of samples}, {starting index}, {ending index} ]
    
    % fill in each row of bFB and uFB with a signal template of the appropriate transit time
    for ii = 1:num_filters
        % get signal templates of the appropriate length
        unipolar_template = imresize( basetemplate_up, [1,transit_times_hr(ii)], 'nearest');
            % row vector of the unipolar signal template [1's and 0's]
        bipolar_template = 2 * unipolar_template - 1;
            % row vector of the bipolar signal template [1's and -1's]
        bw_list(ii) = find(unipolar_template~=unipolar_template(1), 1) - 1;
            % # of [high-resolution] samples in the first segment of the template
        
        % fill the appropriate row of each filterbank by centering the signal template
        start = floor( (maxtime_hr - length(unipolar_template)) / 2 ) + 1; % first index of template
        stop = start + length(unipolar_template) - 1; % last index of template
        uFB_hr(ii,start:stop) = unipolar_template;
        bFB_hr(ii,start:stop) = bipolar_template ./ norm(bipolar_template);
            % scale so that the norm of each template is 1
            % therefore, the correlation won't be larger just because the template length is longer
        
        % fill the appropriate row of the index dictionary
        index_dict_hr(ii,:) = [length(unipolar_template), start, stop];
    end
    
    %% generate low resolution filter bank
    
    % initialize filterbanks, bw_list, & index_dict (low resolution)
    % each row corresponds to a different transit time
    bFB_lr = zeros(num_filters, maxtime_lr); % each row is a template
    uFB_lr = bFB_lr; % each row is a template
    index_dict_lr = zeros(num_filters, 3); % [ {template length in # of samples}, {starting index}, {ending index} ]
    
    for ii = 1:num_filters
        % get signal templates of the appropriate length
        unipolar_template = imresize( basetemplate_up, [1,transit_times_lr(ii)], 'nearest');
            % row vector of the unipolar signal template [1's and 0's]
        bipolar_template = 2 * unipolar_template - 1;
            % row vector of the bipolar signal template [1's and -1's]
        
        % fill the appropriate row of each filterbank by centering the signal template
        start = floor( (maxtime_lr - length(unipolar_template)) / 2 ) + 1; % first index of template
        stop = start + length(unipolar_template) - 1; % last index of template
        uFB_lr(ii,start:stop) = unipolar_template;
        bFB_lr(ii,start:stop) = bipolar_template ./ norm(bipolar_template);
            % scale so that the norm of each template is 1
            % therefore, the correlation won't be larger just because the template length is longer
        
        % fill the appropriate row of the index dictionary
        index_dict_lr(ii,:) = [length(unipolar_template), start, stop];
        
    end
    
    %% get "bandwidth" by inverting segment length [# samples at the resampled rate]
    bw_list = bw_list / (sample_rate/resample_rate);
    bw_list = 1 ./ bw_list;
    
    %% plot, if desired
    if plotting
        plot_FB(bFB_lr, uFB_lr, resample_rate);
        suptitle('low-resolution filter banks');
    end
    
end
