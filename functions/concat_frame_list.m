% concatenate a cell array of detected_events structs from looping over frames with frame_analysis_general
% also concatenates the data and model vectors
% NOTE: frames must be subsequent in the original data, with no gaps or overlaps
% 
% `outfilepath` is an optional string argument (default: empty)
%   > if not empty, results are saved to a file at the provided path
% `plotting` is an optional boolean argument (default: false)
% `D_eff` & `L_ch` & `num_channels` [um,um,int] are optional arguments to calculate cell sizes
%   > if all are provided and `Deff` is not empty, sets `diameter` field in `all_detected_events`
%   > if not, `diameter` is set to `0` for all events

function [all_detected_events, cat_times, cat_data, cat_modelbase, cat_fullmodel, cat_flatmodel, fighandle] = ...
    concat_frame_list(has_cells_list, detected_list, frame_time_list, frame_data_list, ...
                      model_base_list, full_model_list, flat_model_list, ...
                      outfilepath, plotting, D_eff, L_ch, num_channels)

    % input parsing
    if nargin<12
        D_eff = [];

        if nargin<9
            plotting=false;
            if nargin<8
                outfilepath = '';
            end
        end

    elseif isempty(num_channels) || isempty(L_ch)
        D_eff = [];
    end
    
    numframes = length(detected_list);
    %% initialize data vectors
    cat_times = [];
    cat_data = [];
    cat_modelbase = [];
    cat_fullmodel = [];
    cat_flatmodel = [];

    %% loop through frames
    disp('looping through frame results to concatenate...');
    
    for ff = 1:numframes
        
        has_cells = has_cells_list(ff);
        detected_events = detected_list{ff};
        frame_times = frame_time_list{ff};
        frame_data = frame_data_list{ff};
        frame_modelbase = model_base_list{ff};
        frame_fullmodel = full_model_list{ff};
        frame_flatmodel = flat_model_list{ff};

        % get the detected events from the frame
        if has_cells
            if ~exist('all_detected_events','var')
                all_detected_events = detected_events; % this is a deep copy
            else
                all_detected_events = [all_detected_events, detected_events];
            end
        end
        
        % get the data vectors from the frame
        cat_times = [cat_times; frame_times];
        cat_data = [cat_data; frame_data];
        cat_modelbase = [cat_modelbase; frame_modelbase];
        cat_fullmodel = [cat_fullmodel; frame_fullmodel];
        cat_flatmodel = [cat_flatmodel; frame_flatmodel];
        
    end
    
    disp('done looping through frames for concatenation');
    
    %% calculate diameter & rearrange the struct
    disp('finishing concatenated struct...')
    
    if ~exist('all_detected_events','var')
        all_detected_events = [];
        
    else
        
        % Deblois & Bean, if D_eff & L_ch & num_channels were provided
        if isempty(D_eff)
            [all_detected_events.diameter] = deal(0);
        else
            dia = num2cell( calc_diameters(D_eff, all_detected_events, num_channels, L_ch) );
            [all_detected_events.diameter] = dia{:};
            [all_detected_events.D_eff] = deal(D_eff);
            [all_detected_events.L_ch] = deal(L_ch);
            [all_detected_events.num_channels] = deal(num_channels);
        end

        % finalize field order
        if length(fields(all_detected_events)) == 21
            all_detected_events = rmfield(all_detected_events, 'event_num');
            all_detected_events = orderfields(all_detected_events, ...
                {'frame_num', 'frame_event_num', 'sequence_num', 'transit_time', 'amplitude', ...
                'baseline_amp_center', 'start_time','end_time','center_time', ...
                'amp_singlechannel', 'diameter', ...
                'node_amp_factor','node_dur_factor','accel_factor','rise_time', ...
                'detected_centerpos_ix','detected_filterbank_row','detected_correng','detected_pslr','detected_normcorr'});
        elseif length(fields(all_detected_events)) == 24
            all_detected_events = rmfield(all_detected_events, 'event_num');
            all_detected_events = orderfields(all_detected_events, ...
                {'frame_num', 'frame_event_num', 'sequence_num', 'transit_time', 'amplitude', ...
                'baseline_amp_center', 'start_time','end_time','center_time', ...
                'amp_singlechannel', 'diameter', 'D_eff', 'L_ch', 'num_channels', ...
                'node_amp_factor','node_dur_factor','accel_factor','rise_time', ...
                'detected_centerpos_ix','detected_filterbank_row','detected_correng','detected_pslr','detected_normcorr'});
        else
            error('unsupported # of fields in detected event list');
        end
        
    end
    
    disp('done with concatenated event structure');

    %% save concatenated data, if desired
    
    if ~isempty(outfilepath)
        disp('saving results...');
        
        save(outfilepath, 'all_detected_events', 'cat_times', 'cat_data', ...
            'cat_modelbase', 'cat_fullmodel', 'cat_flatmodel');
        
        disp('done saving concatenated results');
    end
    
    %% plot concatenated results, if desired
    if plotting
        disp('plotting...');
        
        fighandle = plot_event_detections(all_detected_events, ...
            cat_times, cat_data, cat_modelbase);
        
        disp('done plotting concatenated results');
    else
        fighandle = [];
    end
    
    %%
    disp('done concatenating frames');
end



