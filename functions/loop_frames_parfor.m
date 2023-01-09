% loop over frames in a parfor using `frame_analysis_general`

%{

INPUTS

* data = (nx1 doubles) = data to be analyzed [units of 1/V or 1/A]
* timedata = (nx1 doubles) = timestamps [sec]
* basedata = (nx1 doubles) = fitted baseline [same units as data]

* noisemask = (nx1 bool) = identifies datapoints that were part of the envelope noise mask in baseline estimation
    > same as noise_mask output from preprocess_data.m
    > used to skip analysis of any frames that contain entirely noise data
    > if you don't have it or don't want to do this check, pass an empty value for this argument

* Nstarts = (nfx1 ints) = start indices for all frames in data
* Nstops = (nfx1 ints) = stop indices for all frames in data
* filterbanks = structure containing the filterbank (& device) information
    * bFB = (mxs doubles) = bipolar filterbank (at data samplerate)
        > same as bFB_ output from generate_FB_general
    * index_dict = (mx3 array) = index information for the filterbanks (at data samplerate)
        > same as index_dict_ output from generate_FB_general
    * transit_times = (mx1 doubles) = list of transit times for the filterbanks [sec]
        > same as transit_times output from generate_FB_general
    * signal_IDs = mx1 column vector specifying which sequence was used for each filterbank row
        > same as signal_IDs output from generate_FB_general
    * signal_code_table = (num_seqsxseq_len matrix of 1's and 0's)
        > each row is one sequence that is encoded in the device geometry (1=pore, 0=node)
        > the row # corresponds to the signal/sequence ID
        > same as `sequences` input to generate_FB_general
    * seqID_perchan = (num_channelsx1 vector of ints)
        > lists the signal/sequence ID used in each channel
        > the vector index corresponds to the channel number

* SIC_param_af = (optional, struct) = parameters for SIC process
    > pass an empty struct to use defaults for all parameters
    > to use non-default values, may include fields:
        * maxSICIter = max # of times to run SIC (default 4)
        * thr = minimum acceptable (multichannel) amplitude of a cell event (default 100)
        * beta = minimum acceptable PSLR for detected event (default 5)
        * norm_corr_thr = minimum acceptable peak normalized correlation for detected event (default 0.2)
* fitting_param_structures = (optional, struct) = initial guesses & parameter bounds for event model fitting
    > pass an empty struct to use all defaults
    > to use non-default values, may include fields:
        * eventfit_initial_guesses = (struct) = initial guesses for event model parameters
            > `initial_guess_struct` input to get_eventlist_params_and_bounds()
        * eventfit_bounds = (struct) = event model parameter bounds, used for final model fit
            > `bounds_struct` input to get_eventlist_params_and_bounds()
        * eventfit_bounds_SIC = (struct) = event model parameter bounds, used for interim model fitting during SIC
            > set to [] to use a copy of eventfit_bounds (default behavior if not provided)
            > set to an empty struct to use all default bounds
        * eventmodel_risetime = (double) [sec] = risetime for piecewise logistic model (default 0.001)
            > time to rise from 10% to 90% of the difference between R_signal_node and R_signal_pore

* SIC_verbose = (optional, bool) = whether to print status updates & results
    > default = false
* SIC_plotting = (optional, bool) = whether to plot the process & results
    > default = false
* use_parallel = (optional, bool)
    > if false, it will not start nor end a parpool and will execute a normal for-loop over the frames
    > if true and there are multiple frames to analyze, it will attempt to 
      start a new parpool, and execute a parfor-loop over the frames
    > default = true
* is_batch = (optional, bool)
    > set to true if function will be called within a batch process
        * in this case, it will not start a new parallel pool
    > if false, it will delete any current parpool and start a new 'threads' parpool
    > default = false

%}

%% output arguments

function [has_cells_list, detected_list, frame_time_list, frame_data_list, model_base_list, full_model_list, flat_model_list] = ...
    loop_frames_parfor(...
        data, timedata, basedata, noisemask, Nstarts, Nstops, filterbanks, ...
        SIC_param_af, fitting_param_structures, ...
        SIC_verbose, SIC_plotting, use_parallel, is_batch)
    
    %% parse inputs

    if nargin<7
        SIC_param_af = struct();
    end
    if nargin<8
        fitting_param_structures = struct();
    end
    if nargin<9
        SIC_verbose = false;
    end
    if nargin<10
        SIC_plotting = false;
    end
    if nargin<11
        use_parallel = true;
    end
    if nargin<12
        is_batch = false;
    end
    
    % check frame breakpoints
    if length(Nstarts) ~= length(Nstops)
        error('Nstarts & Nstops must have the same # of elements')
    end
    
    % don't use parallel if only 1 frame
    numframes = length(Nstarts);
    if numframes==1
        use_parallel = false;
    end

    %% setup

    % extract broadcast variables from structs
    bFB = filterbanks.bFB;
    index_dict = filterbanks.index_dict;
    fb_transittimes = filterbanks.transit_times;
    fb_signal_IDs = filterbanks.signal_IDs;
    signal_code_table = filterbanks.signal_code_table;
    seqID_perchan = filterbanks.seqID_perchan;

    % preallocate lists
    has_cells_list = zeros(1,numframes);
    detected_list = cell(1,numframes);
    frame_time_list = cell(1,numframes);
    frame_data_list = cell(1,numframes);
    model_base_list = cell(1,numframes);
    full_model_list = cell(1,numframes);
    flat_model_list = cell(1,numframes);

    %% loop over the frames
    
    if use_parallel
        %% loop over all the frames in a parfor
        
        if ~is_batch
            delete(gcp('nocreate'));
            parpool('threads');
%             parpool('local', feature('numcores')); % process pool using the # of physical cores
%             parpool('local'); % process pool (might use the # of logical cores instead of the # of physical cores)
        end

        tic;
        parfor ff = 1:numframes
            
            SIC_plotting = false; % assume no SIC plotting when using parfor

            [has_cells, detected_events, frame_time, frame_data, model_base, full_model] = ...
                 frame_analysis_general(...
                    ff, numframes, Nstarts, Nstops, ...
                    data, timedata, basedata, noisemask, signal_code_table, seqID_perchan, ...
                    bFB, index_dict, fb_transittimes, fb_signal_IDs, ...
                    SIC_param_af, fitting_param_structures, SIC_plotting, SIC_verbose);

            % add parsed outputs to overall lists
            has_cells_list(ff) = has_cells;
            detected_list{ff} = detected_events;
            frame_time_list{ff} = frame_time;
            frame_data_list{ff} = frame_data;
            model_base_list{ff} = model_base;
            full_model_list{ff} = full_model;
            flat_model_list{ff} = full_model - model_base;

        end

        fprintf(['******************************** ','finished looping through %u frames!', ...
            ' ********************************\n'], numframes);
        toc;
        
    else
        %% loop over the frames in a normal for-loop
        
        tic;
        for ff = 1:numframes

            [has_cells, detected_events, frame_time, frame_data, model_base, full_model] = ...
                 frame_analysis_general(...
                    ff, numframes, Nstarts, Nstops, ...
                    data, timedata, basedata, noisemask, signal_code_table, seqID_perchan, ...
                    bFB, index_dict, fb_transittimes, fb_signal_IDs, ...
                    SIC_param_af, fitting_param_structures, SIC_plotting, SIC_verbose);

            % add parsed outputs to overall lists
            has_cells_list(ff) = has_cells;
            detected_list{ff} = detected_events;
            frame_time_list{ff} = frame_time;
            frame_data_list{ff} = frame_data;
            model_base_list{ff} = model_base;
            full_model_list{ff} = full_model;
            flat_model_list{ff} = full_model - model_base;

        end

        fprintf(['******************************** ','finished looping through %u frames!', ...
            ' ********************************\n'], numframes);
        toc;
        
    end
    
end









