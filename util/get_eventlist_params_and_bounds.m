% fitting_objective_function takes `params` input in the form of a repeating [6nx1] column vector, where n is the # of events
% lsqnonlin takes bounds inputs `lb` & `ub`, which each should also be [6nx1] vector
% for each event, we need the following parameters in this order:
%   [event_t_centers, event_t_transits, event_amps, event_node_amp_factors, node_dur_factors, accel_factors]
%       > repeat this vector n times, with the params for each event

function [params_vector, lb_vector, ub_vector] = get_eventlist_params_and_bounds(...
    event_list, timestamps, transit_times_fb, initial_guess_struct, bounds_struct)
    
%% deal with input args

    if nargin<4
        initial_guess_struct = [];
    end
    
    if nargin<5
        bounds_struct = [];
    end

    %% defaults: initial guesses

    if isempty(initial_guess_struct)
        initial_guess_struct = struct();
    end

    if ~isfield(initial_guess_struct, 'node_amp_factor')
        % for the data in test_analyze_frame, I found ratios of 0.3445 (smaller event) & 0.2717 (larger event)
        initial_guess_struct.node_amp_factor = 0.3;
    end
    if ~isfield(initial_guess_struct, 'node_dur_factor')
        initial_guess_struct.node_dur_factor = 1.15;
    end
    if ~isfield(initial_guess_struct, 'accel_factor')
        initial_guess_struct.accel_factor = 1;
    end

    %% defaults: bounds

    if isempty(bounds_struct)
        bounds_struct = struct();
    end

    if ~isfield(bounds_struct, 'max_timeshift_sec')
        bounds_struct.max_timeshift_sec = 0.007;
    end
    if ~isfield(bounds_struct, 'max_ttime_scalefactor')
        bounds_struct.max_ttime_scalefactor = 0.1;
    end
    if ~isfield(bounds_struct, 'max_amp_scalefactor')
        bounds_struct.max_amp_scalefactor = 1;
    end
    if ~isfield(bounds_struct, 'max_nodeampfactor_scalefactor')
        bounds_struct.max_nodeampfactor_scalefactor = .5;
    end
    if ~isfield(bounds_struct, 'max_nodedurfactor_scalefactor')
        bounds_struct.max_nodedurfactor_scalefactor = .15;
    end
    if ~isfield(bounds_struct, 'max_accelfactor_scalefactor')
        bounds_struct.max_accelfactor_scalefactor = .3;
    end
    
    %% fill in event parameters

    % 6xn matrix where each column is one event's parameters
    num_events = length(event_list);
    params_matrix = nan(6,num_events);
    for ee = 1:num_events
        params_matrix(1,ee) = timestamps(event_list(ee).pos);
        params_matrix(2,ee) = transit_times_fb(event_list(ee).ttindex);
        params_matrix(3,ee) = event_list(ee).single_channel_amp;
        params_matrix(4,ee) = initial_guess_struct.node_amp_factor;
        params_matrix(5,ee) = initial_guess_struct.node_dur_factor;
        params_matrix(6,ee) = initial_guess_struct.accel_factor;
    end

    % convert the matrix column-major into a vector
    params_vector = params_matrix(:);

    %% fill in parameter bounds per-event

    % 6xn matrices where each column is one event's parameter bounds
    lb_matrix = nan(6,num_events);
    ub_matrix = nan(6,num_events);
    for ee = 1:num_events
        % t_center
        lb_matrix(1,ee) = params_matrix(1,ee) - bounds_struct.max_timeshift_sec;
        ub_matrix(1,ee) = params_matrix(1,ee) + bounds_struct.max_timeshift_sec;

        % transit time
        lb_matrix(2,ee) = params_matrix(2,ee) - params_matrix(2,ee)*bounds_struct.max_ttime_scalefactor;
        ub_matrix(2,ee) = params_matrix(2,ee) + params_matrix(2,ee)*bounds_struct.max_ttime_scalefactor;

        % amplitude
        lb_matrix(3,ee) = params_matrix(3,ee) - params_matrix(3,ee)*bounds_struct.max_amp_scalefactor;
        ub_matrix(3,ee) = params_matrix(3,ee) + params_matrix(3,ee)*bounds_struct.max_amp_scalefactor;
        lb_matrix(3,ee) = max(0, lb_matrix(3,ee)); % ensure amplitude can't be negative

        % node_amp_factor
%         lb_matrix(4,ee) = params_matrix(4,ee) - params_matrix(4,ee)*bounds_struct.max_nodeampfactor_scalefactor;
        ub_matrix(4,ee) = params_matrix(4,ee) + params_matrix(4,ee)*bounds_struct.max_nodeampfactor_scalefactor;
%         lb_matrix(4,ee) = max(0, lb_matrix(4,ee)); % ensure >=0
%         lb_matrix(4,ee) = min(0, lb_matrix(4,ee)); % ensure 0 is included in the range
        lb_matrix(4,ee) = 0; % lower-bound should always be 0
        ub_matrix(4,ee) = min(1, ub_matrix(4,ee)); % ensure <=1

        % node_dur_factor
        lb_matrix(5,ee) = params_matrix(5,ee) - params_matrix(5,ee)*bounds_struct.max_nodedurfactor_scalefactor;
        ub_matrix(5,ee) = params_matrix(5,ee) + params_matrix(5,ee)*bounds_struct.max_nodedurfactor_scalefactor;
        lb_matrix(5,ee) = min(1, lb_matrix(5,ee)); % ensure 1 is included in the range
        ub_matrix(5,ee) = max(1, ub_matrix(5,ee)); % ensure 1 is included in the range

        % accel_factor
        lb_matrix(6,ee) = params_matrix(6,ee) - params_matrix(6,ee)*bounds_struct.max_accelfactor_scalefactor;
        ub_matrix(6,ee) = params_matrix(6,ee) + params_matrix(6,ee)*bounds_struct.max_accelfactor_scalefactor;
        lb_matrix(6,ee) = min(1, lb_matrix(6,ee)); % ensure 1 is included in the range
        ub_matrix(6,ee) = max(1, ub_matrix(6,ee)); % ensure 1 is included in the range
    end

    % convert the matrices column-major into vectors
    lb_vector = lb_matrix(:);
    ub_vector = ub_matrix(:);

end