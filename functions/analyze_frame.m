%{

INPUTS

* frame = (fx1 doubles) = data to be analyzed (units of 1/V or 1/A)
* baseline = (fx1 doubles) = fitted baseline for frame (same units as frame)
* timestamps = (fx1 doubles) = timestamps for frame & baseline, in sec

* seqID_perchan = (num_channelsx1 vector of ints)
    > lists the signal/sequence ID used in each channel
    > the vector index corresponds to the channel number
* signal_code_table = (num_seqsxseq_len matrix of 1's and 0's)
    > each row is one sequence that is encoded in the device geometry (1=pore, 0=node)
    > the row # corresponds to the signal/sequence ID
    > same as `sequences` input to generate_FB_general

* bpFB = (mxs doubles) = bipolar filterbank (at frame samplerate)
    > same as bFB_lr output from generate_FB_general
* index_dict = (mx3 array) = index information for the filterbanks (at frame samplerate)
    > same as index_dict_lr output from generate_FB_general
* signal_IDs = mx1 column vector specifying which sequence was used for each filterbank row
    > same as signal_IDs output from generate_FB_general
* fb_transittimes = (mx1 doubles) = list of transit times for the filterbanks [sec]
    > same as transit_times output from generate_FB_general

* SIC_param = (optional, struct) = parameters for SIC process
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

* plotting = (optional, bool) = whether to plot the process & results
    > default = false
* verbose = (optional, bool) = whether to print status updates & results
    > default = true

OUTPUTS

* detected_events = struct array = final list of detected events
    * event_num = (int) = event ID # within frame
    * amplitude = (double) = event (multichannel) amplitude [same units as frame]
    * baseline_amp_center = (double) = baseline amplitude at center of event [same units as Frame]
    * transit_time = (double) = event transit time [sec]
    * start_time = event start time [sec]
    * end_time = (double) = event end time [sec]
    * center_time = (double) = timestamp at center of event [sec]
    * sequence_num = (int) = event sequence/signal ID
    * amp_singlechannel = (double) = fitted single-channel delta_R
    * node_amp_factor = (double) = fitted ratio (R_signal_node / R_signal_pore)
    * node_dur_factor = (double) = fitted ratio (dT_node / dT_pore)
    * accel_factor = (double) = fitted ratio (V_end / V_start) (how much the cell is speeding up or slowing down during the event)
    * rise_time = (double) = time to rise from 10% to 90% of the difference between R_signal_node and R_signal_pore
    * detected_centerpos_ix = (int) = time-index within frame of the inital matched-filter detection
    * detected_filterbank_row = (int) = row within the filterbank of the inital matched-filter detection
    * detected_correng = (double) = peak correlation energy of the inital matched-filter detection
    * detected_pslr = (double) = calculated PSLR of the inital matched-filter detection
    * detected_normcorr = (double) = normalized peak correlation of the inital matched-filter detection

* final_model = fx1 doubles = fitted model of events (including the baseline) [same units as frame]

NOTES
>> todo: support for sequences having different base times & lengths
>> note: '../Util' must be on search path
    > for apply_FB, adaptiveDetection, and getAll
%}

function [detected_events, final_model] = analyze_frame( ...
    frame, baseline, timestamps, seqID_perchan, signal_code_table, ...
    bpFB, index_dict, signal_IDs, fb_transittimes, ...
    SIC_param, fitting_param_structures, plotting, verbose)

    %% TEMP

    neg_penalty = 2;
    force_non_neg = true;
    improvement_ratio_threshold = 0.05; % proportion of bof_all (should be >0)
    
    %% parse input arguments
    
    p = inputParser;
    p.FunctionName = 'analyze_frame';
    p.PartialMatching = false;
    p.CaseSensitive = true;
    p.StructExpand = false;
    
    % required parameters (data information)
    p.addRequired('frame', @(x) isnumeric(x) && iscolumn(x) && size(x,1)>1 );
    p.addRequired('baseline', @(x) isnumeric(x) && iscolumn(x) && size(x,1)>1 );
    p.addRequired('timestamps', @(x) isnumeric(x) && iscolumn(x) && size(x,1)>1 );
    p.addRequired('seqID_perchan', @(x) isnumeric(x) && (isrow(x) || iscolumn(x)) );
    p.addRequired('signal_code_table', @(x) isnumeric(x) && ismatrix(x) && size(x,1)>=1 && size(x,2)>1 );

    % required parameters (filterbanks)
    p.addRequired('bpFB', @(x) isnumeric(x) && ismatrix(x) && size(x,1)>1 && size(x,2)>1 );
    p.addRequired('index_dict', @(x) isnumeric(x) && ismatrix(x) && size(x,1)>1 && size(x,2)==3 );
    p.addRequired('signal_IDs', @(x) isnumeric(x) && iscolumn(x) && size(x,1)>1 );
    p.addRequired('fb_transittimes', @(x) isnumeric(x) && iscolumn(x) && size(x,1)>1 );
    
    % optional parameters
    p.addOptional('SIC_param', struct());
    p.addOptional('fitting_param_structures', struct());
    p.addOptional('plotting', false, @check_bool);
    p.addOptional('verbose', true, @check_bool);
    
    % run parser based on the provided arguments
    if nargin==9
        parse(p, frame, baseline, timestamps, seqID_perchan, signal_code_table, bpFB, index_dict, signal_IDs, fb_transittimes);
    elseif nargin==10
        parse(p, frame, baseline, timestamps, seqID_perchan, signal_code_table, bpFB, index_dict, signal_IDs, fb_transittimes, SIC_param);
    elseif nargin==11
        parse(p, frame, baseline, timestamps, seqID_perchan, signal_code_table, bpFB, index_dict, signal_IDs, fb_transittimes, SIC_param, fitting_param_structures);
    elseif nargin==12
        parse(p, frame, baseline, timestamps, seqID_perchan, signal_code_table, bpFB, index_dict, signal_IDs, fb_transittimes, SIC_param, fitting_param_structures, plotting);
    elseif nargin >= 13
        if nargin>13
            warning('ignoring extra arguments');
        end
        parse(p, frame, baseline, timestamps, seqID_perchan, signal_code_table, bpFB, index_dict, signal_IDs, fb_transittimes, SIC_param, fitting_param_structures, plotting, verbose);
    else
        error('not enough arguments');
    end
    
    % get flags
    plotting = p.Results.plotting;
    verbose = p.Results.verbose;
    
    % check dimensional agreement between arguments
    if numel(frame) ~= numel(baseline)
        error('frame and baseline must be the same length');
    end
    if numel(frame) ~= numel(timestamps)
        error('frame and timestamps must be the same length');
    end
    if size(bpFB,1) ~= size(index_dict,1)
        error('index_dict must have the same # of rows as filterbank');
    end
    if size(bpFB,1) ~= size(signal_IDs,1)
        error('signal_IDs must have the same # of rows as filterbank');
    end
    if size(bpFB,1) ~= size(fb_transittimes,1)
        error('fb_transittimes must have the same # of rows as filterbank');
    end

    %% parse optional structure arguments

    % get input structs from initial parser (will be empty structs if they weren't provided)
    SIC_param = p.Results.SIC_param;
    fitting_param_structures = p.Results.fitting_param_structures;

    p2 = inputParser;
    p2.FunctionName = 'analyze_frame';
    p2.PartialMatching = false;
    p2.CaseSensitive = true;
    p2.StructExpand = true;

    % optional parameters (from the SIC_param struct)
    p2.addParameter('maxSICIter', 4, @(x) isnumeric(x) && isscalar(x) && x>0 );
    p2.addParameter('thr', 100, @(x) isnumeric(x) && isscalar(x) && x>0 );
    p2.addParameter('beta', 5, @(x) isnumeric(x) && isscalar(x) );
    p2.addParameter('norm_corr_thr', 0.2, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    
    % optional parameters (from fitting_param_structures)
    p2.addParameter('eventfit_initial_guesses', struct());
    p2.addParameter('eventfit_bounds', struct());
    p2.addParameter('eventfit_bounds_SIC', []);
    p2.addParameter('eventmodel_risetime', 0.001, @(x) isnumeric(x) && isscalar(x) && x>0);

    % run parser
    parse(p2, SIC_param, fitting_param_structures);
    
    % get user-defined or default values for SIC parameters from SIC_param
    max_SIC_iter = p2.Results.maxSICIter;
    min_event_amp = p2.Results.thr;
    min_pslr = p2.Results.beta;
    norm_corr_thr = p2.Results.norm_corr_thr;

    % get fitting info from fitting_param_structures
    eventfit_initial_guesses = p2.Results.eventfit_initial_guesses;
    eventfit_bounds = p2.Results.eventfit_bounds;
    eventfit_bounds_SIC = p2.Results.eventfit_bounds_SIC;
    eventmodel_risetime = p2.Results.eventmodel_risetime;

    % set eventfit_bounds_SIC as a copy eventfit_bounds, if it wasn't provided
    if isempty(eventfit_bounds_SIC)
        eventfit_bounds_SIC = eventfit_bounds;
    end
    
    %% setup

    num_chan = length(seqID_perchan); % number of channels in the device
    code_len = size(signal_code_table,2); % length of a manchester-encoded sequence (# of nodes+pores)

    % get/check list of unique signal ID's in the combined filterbank
    unique_seqIDs = unique(signal_IDs, 'stable');
    if max(unique_seqIDs)>length(signal_code_table) || min(unique_seqIDs)<1
        error('valid sequence ID''s must be between 1 and the length of signal_code_table');
    end
    if any(~ismember(seqID_perchan, unique_seqIDs))
        error('each channel''s sequence ID must be found in filterbank signal_IDs');
    end
    
    % get seq info to label graphs
    seq_labels = compose('seq%u', unique_seqIDs);
    fb_seqstart_rows = linspace(1, size(bpFB,1)+1, length(unique_seqIDs)+1);
    fb_seqstart_rows = fb_seqstart_rows(1:end-1);

    % keep original pre-fitted baseline (not re-fit jointly with events)
    model_base = baseline;

    % SIC initialization
    frame_delta = frame-baseline;
%     correlation_maps = []; % stores the (non-normalized) correlation map at each SIC step as pages
%     norm_correlation_maps = []; % stores the normalized correlation map at each SIC step as pages
    detect_SIC_list = []; % event detections for this frame

    %% successive interference cancellation
    for sic_num = 1:max_SIC_iter
        if verbose
            disp(['SIC Iteration ' num2str(sic_num)]);
        end
         %% apply matched filterbank

        % median-filter the residuals (windowsize ~1/2 expected pulse length)
        min_pulse_len = round(min(index_dict(:,1))/code_len);
        frame_delta_med = movmedian(frame_delta, round(min_pulse_len/2));
        if force_non_neg
            frame_delta_med( frame_delta_med<0 ) = 0; % hack to get rid of bad negative residuals
        end

        [norm_corr_map, corr_map] = apply_FB(frame_delta_med, bpFB);
%         correlation_maps(:,:,sic_num) = corr_map;
%         norm_correlation_maps(:,:,sic_num) = norm_corr_map;
        
        if plotting
            f1 = figure;
            ax1 = subplot(5,1,1);
            plot(timestamps, frame_delta, 'k');
            hold(ax1,'on');
            plot(timestamps, frame_delta_med, 'b');
            hold(ax1, 'off');
            xlim([min(timestamps), max(timestamps)]);

            ax2 = subplot(5,1,2);
            imagesc(ax2, timestamps, 1:length(fb_transittimes), corr_map);
            set(ax2,'xlim', get(ax1,'xlim'));
            linkaxes([ax1,ax2], 'x');

            if length(fb_seqstart_rows)>1
                yline(ax2, fb_seqstart_rows(2:end));
                yticks(ax2, fb_seqstart_rows + diff(fb_seqstart_rows(1:2))/2 );
            else
                yticks(ax2, length(fb_transittimes)/2);
            end
            set(ax2.YAxis, TickLabels=seq_labels, TickLength=[0,0]);

            ax3 = subplot(5,1,3);
        end

        %% select the best detection
        % use non-normalized correlation to find largest signal & check PSLR threshold
        % (use normalized correlation to check peak correlation threshold)

        cmap_temp = corr_map; % search space: non-normalized correlation with removed points set to nan
%         poslist = 1:size(cmap_temp,2);
        
        % remove any points from the search space that are already in the frame's detection list
        if ~isempty(detect_SIC_list)
            for ee = 1:length(detect_SIC_list)
                oldevent = detect_SIC_list(ee);
                fb_row = oldevent.ttindex;
                pos_ix = oldevent.pos;

                [~, mainlobe_mask] = find_sidelobe_mask(fb_row,pos_ix, size(cmap_temp), index_dict, signal_IDs, code_len);
                cmap_temp(mainlobe_mask) = nan;

            end
        end

        % remove any points from the search space where norm_corr does not exceed threshold
        % but note that we only want to remove it from the search space, not the PSLR threshold check
        cmap_temp2 = cmap_temp;
        cmap_temp2( norm_corr_map<norm_corr_thr ) = nan;

        % find highest correlation peak with PSLR>threshold(beta) and norm_corr>=norm_corr_thresh
        done_searching = false;
        found_event = false;
        while ~done_searching
            [maxval,lix] = max(cmap_temp2, [], 'all');

            % only consider points where the correlation is >0
            if maxval>0

                [fb_row,pos_ix] = ind2sub(size(cmap_temp2),lix);
                PSLR = compute_pslr2d(fb_row,pos_ix, cmap_temp, index_dict, signal_IDs, code_len);
                peak = maxval;
    
                % check that PSLR exceeds the desired threshold
                done_searching = PSLR >= min_pslr;
                if done_searching
                    found_event = true;
                else
                    % remove this point from the search space
                    cmap_temp(fb_row,pos_ix) = nan;
                    cmap_temp2(fb_row,pos_ix) = nan;
                end

            % none of the remaining points have correlation_energy>0 and norm_corr>norm_corr_thresh
            else
                done_searching = true;
                disp('no more events with corr>0');
            end
        end

        %% add current detection to frame list

        if found_event

            % circle the detection on the correlation heatmap
            if plotting
                hold(ax2,'on');
                plot(ax2, timestamps(pos_ix), fb_row, 'ko', MarkerSize=8, LineWidth=1);
                hold(ax2,'off');
            end

            % check that the normalized correlation exceeds the threshold
            % [note: this should now be taken care of by removing the points that don't exceed this threshold from the correlation search space]
            detection_normcorr = norm_corr_map(fb_row,pos_ix);
            if detection_normcorr >= norm_corr_thr

                % infer initial guess for amplitudes based on correlation value
                ebase = baseline(pos_ix);
                template = bpFB(fb_row,:);
                template_max = max(template);
                num_max_samples = sum(template==template_max);
                mca = peak / template_max / num_max_samples;
                sca = 1 / ( 1/(mca + ebase) - (num_chan-1)/(num_chan*ebase) ) - num_chan*ebase;
    
                currentDetection = struct( ...
                    'pos', pos_ix, ... time position index
                    'ttindex', fb_row, ... template row
                    'amp', mca, ...
                    'correng', peak, ...
                    'pslr', PSLR, ...
                    'normcorr', detection_normcorr, ...
                    'blamp', ebase, ...
                    'seqID', signal_IDs(fb_row), ...
                    'single_channel_amp', sca );
                
                % fill in detectSICList with the identified event
                detect_SIC_list = [detect_SIC_list; currentDetection];
                if verbose
                    disp('Adding current detection to frame event list')
                end

            else % the normalized correlation was below the threshold
                fprintf('stopping bc this detection''s normcorr is below threshold: thr=%f, normcorr=%f\n', norm_corr_thr, detection_normcorr);
                break
            end

        else
            % Stopping Conditions: couldn't find any unique detections
            if verbose
                disp('Cannot find anymore detections');
            end
            break
        end

        % plot current frame's detections
        if plotting            
            figure(f1);
            hold(ax3, 'on');
            for ii=1:length(detect_SIC_list)
                plot(ax3, timestamps(detect_SIC_list(ii).pos), fb_transittimes(detect_SIC_list(ii).ttindex), 'b*', LineWidth=1);
            end
            plot(ax3, timestamps(detect_SIC_list(end).pos), fb_transittimes(detect_SIC_list(end).ttindex), 'ro', MarkerSize=10, LineWidth=1);
            set(ax3, 'xlim', get(ax2,'xlim'));
            linkaxes([ax2,ax3], 'x');
            ylabel('transit time [sec]');
            shg;
        end
        
        %% signal model regression (pulse height & shape parameter estimation)

        % set up parameter initial guesses & bounds for lsqnonlin
        %   > turn detectSICList into the initial guess for the vector (params) that can be passed to fitting_objective_function
        %   > also turn human-readable bounds structure into lb & ub vectors for lsqnonlin
        %   * NOTE: for forward_model & single_event_model, event_amp must be the SINGLE-CHANNEL delta_R
        [event_params_init, lb, ub] = get_eventlist_params_and_bounds(...
            detect_SIC_list, timestamps, fb_transittimes, eventfit_initial_guesses, eventfit_bounds_SIC);

        % perform lsqnonlin
        event_signal_IDs = [detect_SIC_list.seqID];
        fun_l2 = @(x) fitting_objective_function(... % produces residuals; lets lsqnonlin use L2 minimization as usual
            x, frame, timestamps, baseline, seqID_perchan, signal_code_table, event_signal_IDs, eventmodel_risetime);

%         fun_l1 = @(x) sqrt(abs(fun_l2(x))); % hack to force lsqnonlin to use L1 minimization
%         pen_nonneg = @(x) ((x<0).*(-x) + 1) .* x; % hack to penalize negative residuals
        pen_nonneg = @(x) (x<0).*(neg_penalty*x) + (x>=0).*(x); % hack to penalize negative residuals
        fun_l1_nonneg = @(x) sqrt(abs(pen_nonneg(fun_l2(x))));
        current_params_fit = lsqnonlin(fun_l1_nonneg, event_params_init, lb, ub, optimset('display','off'));

        % extract fit results
        current_model = frame - fun_l2(current_params_fit);
        current_params_fit_matrix = reshape(current_params_fit, 6, []);
        % update fitted single-channel amplitudes
        current_event_sca = current_params_fit_matrix(3,:);
        sca_cell = num2cell(current_event_sca);
        [detect_SIC_list.single_channel_amp] = deal(sca_cell{:});
        % infer updated multi-channel delta_R (necessary to check against noise threshold)
        for ii=1:length(detect_SIC_list)
            event = detect_SIC_list(ii);
            bl_singlechan = event.blamp * num_chan;
            R_perchan = bl_singlechan * ones(1, num_chan);
            R_perchan(1) = R_perchan(1) + event.single_channel_amp;
            R_parallel = sum(R_perchan.^-1).^-1;
            detect_SIC_list(ii).amp = R_parallel - event.blamp;
        end

        % plot current iteration's fit to frame
        if plotting
            ax4 = subplot(5,1,4);
            plot(timestamps, frame);
            hold(ax4, 'on');
            plot(timestamps, current_model);
            plot(timestamps, model_base);
            set(ax4,'xlim', get(ax1,'xlim'));
            linkaxes([ax1,ax4], 'x');
            shg;
        end

        %% interference cancellation (subtract current model of detected events)

        frame_delta = frame - current_model;

        if plotting
            figure(f1);
            ax5 = subplot(5,1,5);
            plot(timestamps, frame_delta);
            set(ax5,'xlim', get(ax1,'xlim'));
            linkaxes([ax1,ax5], 'x');
            title('residuals');
            shg;
        end

        % stop SIC if any of the fit detections are too small
        excluded_ix = [detect_SIC_list.amp] < min_event_amp;
        if any(excluded_ix)
            excluded_amps = [detect_SIC_list(excluded_ix).amp];
            fprintf('removing detection(s) & stopping b/c below noise floor: thr=%g, x=%s\n', min_event_amp, mat2str(excluded_amps));
            break
        end

    end

    % remove events whose amplitude is below the threshold
    if ~isempty(detect_SIC_list)
        pruned_list_og = detect_SIC_list(~excluded_ix);
    else
        pruned_list_og = [];
    end
    
    %% compute joint model fitting with all SIC event detections
    % make sure each event is actually important for the fit

        done_fitting = false;
        pruned_list = pruned_list_og;

        while ~done_fitting

            if length(pruned_list) <= 1
                done_fitting = true;
            else
        
                %% fit model with all the events
    
                % turn pruned_list into the initial guess for the vector (params) that can be passed to fitting_objective_function
                % also turn human-readable bounds structure into lb & ub vectors for lsqnonlin
                % NOTE: for forward_model & single_event_model, event_amp must be the SINGLE-CHANNEL delta_R
                ebase = baseline([pruned_list.pos])';
                sca = 1 ./ ( 1./([pruned_list.amp] + ebase) - (num_chan-1) ./ (num_chan .* ebase) ) - num_chan .* ebase;
                sca_cell = num2cell(sca);
                [pruned_list.single_channel_amp] = deal(sca_cell{:});
                [event_params_init, lb, ub] = get_eventlist_params_and_bounds(...
                    pruned_list, timestamps, fb_transittimes, eventfit_initial_guesses, eventfit_bounds);
        
                % perform lsqnonlin
                event_signal_IDs = [pruned_list.seqID];
                fun_l2 = @(x) fitting_objective_function(... % produces residuals; lets lsqnonlin use L2 minimization as usual
                    x, frame, timestamps, baseline, seqID_perchan, signal_code_table, event_signal_IDs, eventmodel_risetime);
                fun_l1 = @(x) sqrt(abs(fun_l2(x))); % hack to force lsqnonlin to use L1 minimization
                [event_params_fit, bof_all] = lsqnonlin(fun_l1, event_params_init, lb, ub, optimset('display','off'));
    
                %% fit model, dropping each event
                
                bof_drop_list = nan(1,length(pruned_list));
                for ee=1:length(pruned_list)
                    pl_drop = pruned_list;
                    pl_drop(ee) = [];
        
                    ebase = baseline([pl_drop.pos])';
                    sca = 1 ./ ( 1./([pl_drop.amp] + ebase) - (num_chan-1) ./ (num_chan .* ebase) ) - num_chan .* ebase;
                    sca_cell = num2cell(sca);
                    [pl_drop.single_channel_amp] = deal(sca_cell{:});
                    [event_params_init, lb, ub] = get_eventlist_params_and_bounds(...
                        pl_drop, timestamps, fb_transittimes, eventfit_initial_guesses, eventfit_bounds);
            
                    % perform lsqnonlin
                    event_signal_IDs = [pl_drop.seqID];
                    fun_l2 = @(x) fitting_objective_function(... % produces residuals; lets lsqnonlin use L2 minimization as usual
                        x, frame, timestamps, baseline, seqID_perchan, signal_code_table, event_signal_IDs, eventmodel_risetime);
                    fun_l1 = @(x) sqrt(abs(fun_l2(x))); % hack to force lsqnonlin to use L1 minimization
                    [event_params_fit, bof_drop] = lsqnonlin(fun_l1, event_params_init, lb, ub, optimset('display','off'));
    
                    bof_drop_list(ee) = bof_drop;
                end
    
                %% see if we need to drop an event
    
                improvement_when_dropped = (bof_all - bof_drop_list)/bof_all;
    
                if any(improvement_when_dropped > -1*improvement_ratio_threshold)
                    [~, most_improved_ee] = max(improvement_when_dropped);
                    
                    disp('dropped an event bc it didn''t improve the model fit:');
                    disp(pruned_list(most_improved_ee));
                    fprintf('* bof_drop = %f, bof_all = %f, improvement_ratio = %f\n', ...
                        bof_drop_list(most_improved_ee), bof_all, improvement_when_dropped(most_improved_ee) );

                    pruned_list(most_improved_ee) = [];
                else
                    done_fitting = true;
                end

            end

        end

    %% get super-final model fit, minus any events that didn't improve the fit
    % also, if needed, remove any events with refitted amplitude below the noise threshold & fit again

    if ~isempty(pruned_list)

        done_fitting = false;
        while ~done_fitting
            if isempty(pruned_list)
                disp('no events remaining in pruned_list');
                done_fitting = true;
            else
                % turn pruned_list into the initial guess for the vector (params) that can be passed to fitting_objective_function
                % also turn human-readable bounds structure into lb & ub vectors for lsqnonlin
                % NOTE: for forward_model & single_event_model, event_amp must be the SINGLE-CHANNEL delta_R
                ebase = baseline([pruned_list.pos])';
                sca = 1 ./ ( 1./([pruned_list.amp] + ebase) - (num_chan-1) ./ (num_chan .* ebase) ) - num_chan .* ebase;
                sca_cell = num2cell(sca);
                [pruned_list.single_channel_amp] = deal(sca_cell{:});
                [event_params_init, lb, ub] = get_eventlist_params_and_bounds(...
                    pruned_list, timestamps, fb_transittimes, eventfit_initial_guesses, eventfit_bounds);
        
                % perform lsqnonlin
                event_signal_IDs = [pruned_list.seqID];
                fun_l2 = @(x) fitting_objective_function(... % produces residuals; lets lsqnonlin use L2 minimization as usual
                    x, frame, timestamps, baseline, seqID_perchan, signal_code_table, event_signal_IDs, eventmodel_risetime);
                fun_l1 = @(x) sqrt(abs(fun_l2(x))); % hack to force lsqnonlin to use L1 minimization
                event_params_fit = lsqnonlin(fun_l1, event_params_init, lb, ub, optimset('display','off'));

                % check if any of the newly-refitted amps are below the noise threshold
                event_params_fit_matrix = reshape(event_params_fit, 6, []);
                badix = zeros(1, size(event_params_fit_matrix,2));
                badamps = [];
                for ii=1:size(event_params_fit_matrix,2)
                    event = event_params_fit_matrix(:,ii);

                    [~,centerpos_ix] = min(abs( timestamps-event(1) )); % event(1) = centertime [sec]
                        % # of samples from the start of the filtered data window
                    blamp = model_base(centerpos_ix); % baseline amp at the center position
                    
                    % turn the single-channel amplitude into the multi-channel delta_R
                    R_perchan = blamp*num_chan * ones(1, num_chan); % blamp*num_chan = bl_singlechan
                    R_perchan(1) = R_perchan(1) + event(3); % event(3) = deltaR_singlechan (single-channel "event amplitude")
                    deltaR_multichan = sum(R_perchan.^-1).^-1 - blamp; % sum(R_perchan.^-1).^-1 = R_parallel

                    if deltaR_multichan < min_event_amp
                        badix(ii) = 1;
                        bad_amps = [badamps, deltaR_multichan];
                    end
                end
    
                % remove any re-fitted detections that are too small
                if any(badix)
                    fprintf('removing re-fitted detection(s) that are below noise floor: thr=%g, x=%s\n', min_event_amp, mat2str(bad_amps));
                    pruned_list = pruned_list(~badix);
                else
                    done_fitting = true;
                end

            end
        end

    
    else
        % pruned_list was already empty, before doing the super-final fit
        event_params_fit = [];
        disp('no events in pruned_list')
    end

    % get the final model
    if ~isempty(event_params_fit)
        final_model = frame - fun_l2(event_params_fit);
        event_params_fit_matrix = reshape(event_params_fit, 6, []);
    else
        final_model = model_base;
    end

    % plot the final model
    if plotting && ~isempty(event_params_fit)
        figure;

        ax1 = subplot(2,1,1);
        plot(frame);
        hold(ax1, 'on');
        plot(final_model);
        title('Model Fit');

        ax2 = subplot(2,1,2);
        plot(frame - final_model);
        title('Residual');

        linkaxes([ax1,ax2], 'x');
    end

    % print final fitting results
    if verbose
        if isempty(event_params_fit)
            disp('detected no particles');
        else
            disp('Final Detected Particle Pulse Heights (single-channel amplitudes): ');
            disp(num2str( event_params_fit_matrix(3,:) ));
        end
    end
    
    %% parse out the final event detections for this frame
    % go through each event in event_params_fit and get the info
    % save important properties to a new struct array

    detected_events = struct();
    
    if ~isempty(event_params_fit)
        for ii=1:size(event_params_fit_matrix,2)

            event = event_params_fit_matrix(:,ii);

            deltaR_singlechan = event(3); % single-channel "event amplitude"
            transittime = event(2); % [sec]

            % figure out the baseline at the event center (based on nearest timestamp index)
            centertime = event(1); % [sec]
            [~,centerpos_ix] = min(abs( timestamps-centertime ));
                % # of samples from the start of the filtered data window
            blamp = model_base(centerpos_ix); % baseline amp at the center position

            % turn the single-channel amplitude into the multi-channel delta_R
            bl_singlechan = blamp * num_chan;
            R_perchan = bl_singlechan * ones(1, num_chan);
            R_perchan(1) = R_perchan(1) + deltaR_singlechan;
            R_parallel = sum(R_perchan.^-1).^-1;
            deltaR_multichan = R_parallel - blamp;

            % fill in detected_events structure

            detected_events(ii).event_num = ii;
            detected_events(ii).amplitude = deltaR_multichan;
            detected_events(ii).baseline_amp_center = blamp;
            detected_events(ii).transit_time = transittime; % [sec]
            detected_events(ii).start_time = centertime - transittime/2; % [sec]
            detected_events(ii).end_time = centertime + transittime/2; % [sec]
            detected_events(ii).center_time = centertime; % [sec]
            
            detected_events(ii).sequence_num = pruned_list(ii).seqID;

            detected_events(ii).amp_singlechannel = deltaR_singlechan;
            detected_events(ii).node_amp_factor = event(4);
            detected_events(ii).node_dur_factor = event(5);
            detected_events(ii).accel_factor = event(6);
            detected_events(ii).rise_time = eventmodel_risetime;

            % fill in info about initial matched-filter detection
            detected_events(ii).detected_centerpos_ix = pruned_list(ii).pos; % time-index within frame of initial detection
            detected_events(ii).detected_filterbank_row = pruned_list(ii).ttindex; % filterbank row of initial detection
            detected_events(ii).detected_correng = pruned_list(ii).correng; % [non-normalized] peak correlation energy of initial detection
            detected_events(ii).detected_pslr = pruned_list(ii).pslr; % PSLR of initial detection [from non-normalized map]
            detected_events(ii).detected_normcorr = pruned_list(ii).normcorr; % normalized peak correlation energy of initial detection

        end
    end
    
    % plot the events detected in this frame (with event #'s, without saving)
    if plotting
        plot_event_detections(detected_events, timestamps, frame, model_base, true);
    end
    
end