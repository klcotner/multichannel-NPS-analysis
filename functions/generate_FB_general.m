% make a filter bank including templates from one or more codes

%{
INPUTS
* sequences = (num_seqsxseq_len matrix of 1's and 0's)
    > each row is one sequence that is encoded in the device geometry (1=pore, 0=node)
* sample_rate = acquisition sample rate, in Hz
* resample_rate = desired sample rate after resampling, in Hz
* FB_param = (optional, struct) = filter bank information
    > uses same FB_param for all sequences
    > to use non-default values, may include fields:
        * num_filters = number of signal templates to generate in the filterbank (default = 100)
        * base_time = expected transit time of one event [sec] (default = 0.2 sec)
        * min_ratio = shortest signal template to generate, as a fraction of base_time (default = 0.05)
        * max_ratio = longest signal template to generate, as a fraction of base_time (default = 1.5)
	> if FB_param is not provided, default values are used for all 4 parameters
* plotting = (optional, bool) = whether or not to plot the combined filterbank
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
* signal_IDs = mx1 column vector specifying which sequence was used for the corresponding filterbank rows

NOTES
>> todo: if max_ratio*base_time is not the same, FB's may have templates of different lengths
>> todo: support for calibrated rcode
>> note: '../Util' must be on the search path
    > for generate_FB & plot_FB_general
%}

function [bFB_hr, uFB_hr, index_dict_hr, bFB_lr, uFB_lr, index_dict_lr, transit_times, bw_list, signal_IDs] = ...
    generate_FB_general(sequences, sample_rate, resample_rate, FB_param, plotting)

    %% parse inputs
    if nargin<5
        plotting = false;
        
        if nargin<4
            FB_param = {};
            
        end
    end
    
    num_seqs = size(sequences,1);
    seq_len = size(sequences,2);
    
    %% preallocate output lists
    % m = maximum # of samples
    % n = # of templates per sequence
    
    % high-res nxm arrays 
    bFB_hr_list = cell(num_seqs,1);
    uFB_hr_list = cell(num_seqs,1);
    
    % high-res nx3 arrays
    index_dict_hr_list = cell(num_seqs,1);
    
    % low-res nxm arrays
    bFB_lr_list = cell(num_seqs,1);
    uFB_lr_list = cell(num_seqs,1);
    
    % low-res nx3 arrays
    index_dict_lr_list = cell(num_seqs,1);
    
    % nx1 column vectors
    transit_times_list = cell(num_seqs,1);
    bw_list_list = cell(num_seqs,1);
    signal_IDs_list = cell(num_seqs,1);
    
    
    %% generate each filter bank (don't plot yet)
    
    for seqix = 1:num_seqs
    
        % generate single-sequence filterbank
        [bFB_hr_i, uFB_hr_i, index_dict_hr_i, bFB_lr_i, uFB_lr_i, index_dict_lr_i, transit_times_i, bw_list_i] = ...
            generate_FB(sequences(seqix,:), sample_rate, resample_rate, FB_param, false);
        
        % set the signal IDs for this sequence
        signal_IDs_i = seqix * ones(length(transit_times_i), 1);
        
        % add the single-sequence outputs to overall lists
        bFB_hr_list{seqix} = bFB_hr_i;
        uFB_hr_list{seqix} = uFB_hr_i;
        index_dict_hr_list{seqix} = index_dict_hr_i;
        bFB_lr_list{seqix} = bFB_lr_i;
        uFB_lr_list{seqix} = uFB_lr_i;
        index_dict_lr_list{seqix} = index_dict_lr_i;
        transit_times_list{seqix} = transit_times_i;
        bw_list_list{seqix} = bw_list_i;
        signal_IDs_list{seqix} = signal_IDs_i;
        
    end
    
    %% combine the filter banks
    
    % high-res nxm arrays
    bFB_hr = cell2mat(bFB_hr_list);
    uFB_hr = cell2mat(uFB_hr_list);
    index_dict_hr = cell2mat(index_dict_hr_list);
    
    % low-res nxm arrays
    bFB_lr = cell2mat(bFB_lr_list);
    uFB_lr = cell2mat(uFB_lr_list);
    index_dict_lr = cell2mat(index_dict_lr_list);
    
    % nx1 column vectors
    transit_times = cell2mat(transit_times_list);
    bw_list = cell2mat(bw_list_list);
    signal_IDs = cell2mat(signal_IDs_list);
    
    
    %% plot, if desired
    
    if plotting
        
        % infer the number of templates per sequence (they all use the same FB_param)
        num_templates = size(transit_times_list{1}, 1);
        seq_row_starts = (0:num_seqs-1) * num_templates + 1;
        
        plot_FB_general(bFB_lr, uFB_lr, resample_rate, seq_row_starts);
        
        sgtitle('low-resolution filter banks');
    end
    
end
