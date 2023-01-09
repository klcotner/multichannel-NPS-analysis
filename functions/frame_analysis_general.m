% this function can be called inside of a parfor loop across frames

%{
`plotting` and `verbose` are optional boolean arguments (default: false)

INPUTS
* ff = (int) = index of the frame to analyze
* numframes = (int) = total number of frames that will be analyzed in the enclosing loop (used for progress display)
* Nstarts = (nfx1 ints) = start indices for all frames in data
* Nstops = (nfx1 ints) = stop indices for all frames in data

* data = (nx1 doubles) = data to be analyzed [units of 1/V or 1/A]
* timedata = (nx1 doubles) = timestamps [sec]
* basedata = (nx1 doubles) = fitted baseline [same units as data]

* noisemask = (nx1 bool) = identifies datapoints that were part of the envelope noise mask in baseline estimation
    > same as noise_mask output from preprocess_data.m
    > used to skip analysis of any frames that contain entirely noise data
    > if you don't have it or don't want to do this check, pass an empty value for this argument

* signal_code_table = (num_seqsxseq_len matrix of 1's and 0's)
    > each row is one sequence that is encoded in the device geometry (1=pore, 0=node)
    > the row # corresponds to the signal/sequence ID
    > same as `sequences` input to generate_FB_general
* seqID_perchan = (num_channelsx1 vector of ints)
    > lists the signal/sequence ID used in each channel
    > the vector index corresponds to the channel number

* bFB = (mxs doubles) = bipolar filterbank (at data samplerate)
    > same as bFB_ output from generate_FB_general
* index_dict = (mx3 array) = index information for the filterbanks (at data samplerate)
    > same as index_dict_ output from generate_FB_general
* fb_transittimes = (mx1 doubles) = list of transit times for the filterbanks [sec]
    > same as transit_times output from generate_FB_general
* fb_signal_IDs = mx1 column vector specifying which sequence was used for each filterbank row
    > same as signal_IDs output from generate_FB_general

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

* SIC_plotting = (optional, bool) = whether to plot the process & results
    > default = false
* SIC_verbose = (optional, bool) = whether to print status updates & results
    > default = true
    
%}

function [has_cells, detected_events, frame_time, frame_data, model_base, full_model] = ...
         frame_analysis_general(...
            ff, numframes, Nstarts, Nstops, ...
            data, timedata, basedata, noisemask, signal_code_table, seqID_perchan, ...
            bFB, index_dict, fb_transittimes, fb_signal_IDs, ...
            SIC_param_af, fitting_param_structures, SIC_plotting, SIC_verbose)

        % default parameter options
        if nargin<14
            SIC_param_af = struct();
        end
        if nargin<15
            fitting_param_structures = struct();
        end

        % plotting & verbose are false by default
        if nargin<16
            SIC_plotting = false;
        end
        if nargin<17
            SIC_verbose = false;
        end

        fprintf(['******************************** FRAME # %u / %u ',...
            '********************************\n'], ff, numframes);

        % check if frame is noise-only
        ixlist = Nstarts(ff) : Nstops(ff);
        if isempty(noisemask) % if noisemask argument wasn't provided, then skip this check
            skip_frame_detection = false;
        else
            skip_frame_detection = all(noisemask(ixlist));
        end

        % extract frame data
        frame_data = data(ixlist);
        frame_time = timedata(ixlist);
        model_base = basedata(ixlist); % don't redo ASLS on each frame (baseline doesn't get re-fitted with each frame)

        if ~skip_frame_detection
            % do SIC using analyze_frame
            [detected_events, full_model] = analyze_frame( ...
                frame_data, model_base, frame_time, seqID_perchan, signal_code_table, ...
                bFB, index_dict, fb_signal_IDs, fb_transittimes, ...
                SIC_param_af, fitting_param_structures, SIC_plotting, SIC_verbose);
    
            % add frame information to detected_events
            has_cells = ~isempty(detected_events) && ~isempty(struct2cell(detected_events));
            if has_cells
                [detected_events.frame_num] = deal(ff);
                [detected_events.frame_event_num] = deal(detected_events.event_num);
            end

        % if skipping this frame, assume no cells & fill in the model with the original baseline
        else
            has_cells = false;
            detected_events = [];
            full_model = model_base;
        end
    
end