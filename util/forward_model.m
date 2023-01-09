% forward model for parallel channels, with shape parameters for the signal model
% 
% INPUTS:
%   * timestamps - frame timestamps [sec] (frame_len x 1)
%   * baseline - full baseline resistance [Ohm] (frame_len x 1)
%       > NOTE: assumes each channel has the same baseline resistance equal to total baseline x num_chans
%   * signal_ID_per_chan = which signal_ID goes with which channel (num_channels x 1)
%   * signal_code_table = [num_codes x code_length] each row contains the code for the corresponding signal_ID
% ...
%   * event_signal_IDs = vector of signal_ID for each detected event
%   * event_t_centers = vector of t_center for each detected event [sec]
%   * event_t_transits = vector of t_transit for each detected event [sec]
%   * event_amps = vector of amplitudes for each detected event (single-channel delta_R in the pores) [ohms]
%       ** this should be delta_R WITHIN JUST ONE CHANNEL
%   * event_node_amp_factors = vector of node_amp_factor for each detected event (single-channel ratio of delta_R_node/delta_R_pore)
%   * event_node_dur_factors = vector of node_dur_factor for each detected event (ratio of dT_node/dT_pore)
%   * event_accel_factors = vector of accel_factor for each detected event (ratio V_end/V_start)
%
% OUTPUT:
%   * y_hat - estimated total (parallel) resistance [Ohm]

function y_hat = forward_model(timestamps, baseline, signal_ID_per_chan, signal_code_table, ...
    event_signal_IDs, event_t_centers, event_t_transits, event_amps, ...
    event_node_amp_factors, event_node_dur_factors, event_accel_factors, rise_time)
    
    %% setup
    num_chans = length(signal_ID_per_chan); % infer number of channels
    num_events = length(event_signal_IDs);
    frame_len = length(timestamps);
    
    b_chan = baseline * num_chans; % assumes each channel's baseline resistance is equal to total baseline resistance x num_chans
    y_hat_per_chan = zeros(frame_len, num_chans); % estimated resistance per channel vs. time
    for i = 1:num_chans
        y_hat_per_chan(:,i) = b_chan;
    end
    
    %% assign events to channels
    % NOTE: assumes all events of the same signal_ID came from the same channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO: reconsider this assumption
    chan_per_event = nan(size(event_signal_IDs));
    for i = 1:num_events
        matching_chan = find(signal_ID_per_chan == event_signal_IDs(i), 1, 'first'); % finds first channel index that matched signal_ID of ith event
        chan_per_event(i) = matching_chan;
    end
    
    %% iterate through each event and add it to y_hat_per_chan
    for event_ix = 1:num_events
        chan_num = chan_per_event(event_ix);
    
        event_model = single_event_model(signal_code_table(event_signal_IDs(event_ix),:), ...
            timestamps, ...
            event_t_centers(event_ix), ...
            event_t_transits(event_ix), ...
            event_amps(event_ix), ...
            event_node_amp_factors(event_ix), ...
            event_node_dur_factors(event_ix), ...
            event_accel_factors(event_ix), ...
            rise_time);
    
        y_hat_per_chan(:,chan_num) = y_hat_per_chan(:,chan_num) + event_model;
    end
    
    %% parallel combine
    y_hat = sum(y_hat_per_chan.^-1, 2).^-1;

end
